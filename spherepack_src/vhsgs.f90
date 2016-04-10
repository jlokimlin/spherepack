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
! ... file vhsgs.f
!
!     this file contains code and documentation for subroutines
!     vhsgs and vhsgsi
!
! ... files which must be loaded with vhsgs.f
!
!     sphcom.f, hrfft.f, gaqd.f
!
!     subroutine vhsgs(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, 
!    +                 mdab, ndab, wvhsgs, lvhsgs, work, lwork, ierror)
!                                                                              
!   
!     subroutine vhsgs performs the vector spherical harmonic synthesis
!     of the arrays br, bi, cr, and ci and stores the result in the
!     arrays v and w.  the synthesis is performed on an equally spaced
!     longitude grid and a gaussian colatitude grid (measured from
!     the north pole). v(i, j) and w(i, j) are the colatitudinal and
!     east longitudinal components respectively, located at the i(th)
!     colatitude gaussian point (see nlat below) and longitude
!     phi(j) = (j-1)*2*pi/nlon.  the spectral respresentation of (v, w)
!     is given below at output parameters v, w.
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
!     nt     the number of syntheses.  in the program that calls vhsgs, 
!            the arrays v, w, br, bi, cr, and ci can be three dimensional
!            in which case multiple syntheses will be performed.
!            the third index is the synthesis index which assumes the 
!            values k=1, ..., nt.  for a single synthesis set nt=1. the
!            discription of the remaining parameters is simplified
!            by assuming that nt=1 or that all the arrays are two
!            dimensional.
!
!     idvw   the first dimension of the arrays v, w as it appears in
!            the program that calls vhags. if ityp .le. 2 then idvw
!            must be at least nlat.  if ityp .gt. 2 and nlat is
!            even then idvw must be at least nlat/2. if ityp .gt. 2
!            and nlat is odd then idvw must be at least (nlat+1)/2.
!
!     jdvw   the second dimension of the arrays v, w as it appears in
!            the program that calls vhsgs. jdvw must be at least nlon.
!
!     br, bi  two or three dimensional arrays (see input parameter nt)
!     cr, ci  that contain the vector spherical harmonic coefficients
!            in the spectral representation of v(i, j) and w(i, j) given
!            below at the discription of output parameters v and w.
!
!     mdab   the first dimension of the arrays br, bi, cr, and ci as it
!            appears in the program that calls vhsgs. mdab must be at
!            least min(nlat, nlon/2) if nlon is even or at least
!            min(nlat, (nlon+1)/2) if nlon is odd.
!
!     ndab   the second dimension of the arrays br, bi, cr, and ci as it
!            appears in the program that calls vhsgs. ndab must be at
!            least nlat.
!
!     wvhsgs an array which must be initialized by subroutine vhsgsi.
!            once initialized, wvhsgs can be used repeatedly by vhsgs
!            as long as nlon and nlat remain unchanged.  wvhsgs must
!            not be altered between calls of vhsgs.
!
!     lvhsgs the dimension of the array wvhsgs as it appears in the
!            program that calls vhsgs. define
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
!                 l1*l2*(nlat+nlat-l1+1)+nlon+15+2*nlat
!
!
!     work   a work array that does not have to be saved.
!
!     lwork  the dimension of the array work as it appears in the
!            program that calls vhsgs. define
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
!            v(i, j), w(i, j) contain the components at the guassian colatitude
!            point theta(i) and longitude phi(j) = (j-1)*2*pi/nlon.
!            the index ranges are defined above at the input parameter
!            ityp. v and w are computed from the formulas given below.
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
!     7.   theta(i) = i(th) gaussian grid point and phi(j) = (j-1)*2*pi/nlon
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
!            = 9  error in the specification of lvhsgs
!            = 10 error in the specification of lwork
!
!
!     subroutine vhsgsi(nlat, nlon, wvhsgs, lvhsgs, dwork, ldwork, ierror)
!
!     subroutine vhsgsi initializes the array wvhsgs which can then be
!     used repeatedly by subroutine vhsgs until nlat or nlon is changed.
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
!            nlon = 72 for a five degree grid. nlon must be greater
!            than zero. the axisymmetric case corresponds to nlon=1.
!            the efficiency of the computation is improved when nlon
!            is a product of small prime numbers.
!
!     lvhsgs the dimension of the array wvhsgs as it appears in the
!            program that calls vhsgs. define
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
!                 l1*l2*(nlat+nlat-l1+1)+nlon+15+2*nlat
!
!     dwork a real work array that does not need to be saved
!
!     ldwork the dimension of the array dwork as it appears in the
!            program that calls vhsgsi. ldwork must be at least
!
!                 (3*nlat*(nlat+3)+2)/2

!
!     **************************************************************
!
!     output parameters
!
!     wvhsgs an array which is initialized for use by subroutine vhsgs.
!            once initialized, wvhsgs can be used repeatedly by vhsgs
!            as long as nlat and nlon remain unchanged.  wvhsgs must not
!            be altered between calls of vhsgs.
!
!
!     ierror = 0  no errors
!            = 1  error in the specification of nlat
!            = 2  error in the specification of nlon
!            = 3  error in the specification of lvhsgs
!            = 4  error in the specification of lwork
!
subroutine vhsgs(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, &
    mdab, ndab, wvhsgs, lvhsgs, work, lwork, ierror)
    dimension v(idvw, jdvw, 1), w(idvw, jdvw, 1), br(mdab, ndab, 1), &
        bi(mdab, ndab, 1), cr(mdab, ndab, 1), ci(mdab, ndab, 1), &
        work(1), wvhsgs(1)
    ierror = 1
    if (nlat < 3) return
    ierror = 2
    if (nlon < 1) return
    ierror = 3
    if (ityp<0 .or. ityp>8) return
    ierror = 4
    if (nt < 0) return
    ierror = 5
    imid = (nlat+1)/2
    if ((ityp<=2 .and. idvw<nlat) .or. &
        (ityp>2 .and. idvw<imid)) return
    ierror = 6
    if (jdvw < nlon) return
    ierror = 7
    mmax = min(nlat, (nlon+1)/2)
    if (mdab < mmax) return
    ierror = 8
    if (ndab < nlat) return
    ierror = 9
    idz = (mmax*(nlat+nlat-mmax+1))/2
    lzimn = idz*imid
    if (lvhsgs < lzimn+lzimn+nlon+15) return
    ierror = 10
    idv = nlat
    if (ityp > 2) idv = imid
    lnl = nt*idv*nlon
    if (lwork < lnl+lnl+idv*nlon) return
    ierror = 0
    ist = 0
    if (ityp <= 2) ist = imid
    !
    !     set wvhsgs pointers
    !
    lmn = nlat*(nlat+1)/2
    jw1 = 1
    jw2 = jw1+imid*lmn
    jw3 = jw2+imid*lmn
    !
    !     set work pointers
    !
    iw1 = ist+1
    iw2 = lnl+1
    iw3 = iw2+ist
    iw4 = iw2+lnl

    call vhsgs1(nlat, nlon, ityp, nt, imid, idvw, jdvw, v, w, mdab, ndab, &
        br, bi, cr, ci, idv, work, work(iw1), work(iw2), work(iw3), &
        work(iw4), idz, wvhsgs(jw1), wvhsgs(jw2), wvhsgs(jw3))

end subroutine vhsgs



subroutine vhsgs1(nlat, nlon, ityp, nt, imid, idvw, jdvw, v, w, mdab, &
    ndab, br, bi, cr, ci, idv, ve, vo, we, wo, work, idz, vb, wb, wrfft)
    dimension v(idvw, jdvw, 1), w(idvw, jdvw, 1), br(mdab, ndab, 1), &
        bi(mdab, ndab, 1), cr(mdab, ndab, 1), ci(mdab, ndab, 1), &
        ve(idv, nlon, 1), vo(idv, nlon, 1), we(idv, nlon, 1), &
        wo(idv, nlon, 1), work(1), wrfft(1), &
        vb(imid, 1), wb(imid, 1)
    nlp1 = nlat+1
    mlat = mod(nlat, 2)
    mlon = mod(nlon, 2)
    mmax = min(nlat, (nlon+1)/2)
    imm1 = imid
    if (mlat /= 0) then
        imm1 = imid-1
    end if

    do k=1, nt
        do j=1, nlon
            do i=1, idv
                ve(i, j, k) = 0.0
                we(i, j, k) = 0.0
            end do
        end do
    end do

    ndo1 = nlat
    ndo2 = nlat

    if (mlat /= 0) then
        ndo1 = nlat-1
    end if

    if (mlat == 0) then
        ndo2 = nlat-1
    end if

    vector_symmetries: select case (ityp)
        case (0)
            !
            !==> case ityp=0   no symmetries
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
            !==> case m = 1 through nlat-1
            !
            if (mmax < 2) then
                exit vector_symmetries
            end if

            do mp1=2, mmax
                m = mp1-1
                !     mb = m*(nlat-1)-(m*(m-1))/2
                mb = m*nlat-(m*(m+1))/2
                mp2 = mp1+1

                if (mp1 > ndo1) then
                    go to 26
                end if

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

                        if (mlat == 0) then
                            exit
                        end if

                        ve(imid, 2*mp1-2, k) = ve(imid, 2*mp1-2, k) &
                            -ci(mp1, np1, k)*wb(imid, mn)
                        ve(imid, 2*mp1-1, k) = ve(imid, 2*mp1-1, k) &
                            +cr(mp1, np1, k)*wb(imid, mn)
                        we(imid, 2*mp1-2, k) = we(imid, 2*mp1-2, k) &
                            -bi(mp1, np1, k)*wb(imid, mn)
                        we(imid, 2*mp1-1, k) = we(imid, 2*mp1-1, k) &
                            +br(mp1, np1, k)*wb(imid, mn)
                    end do
                end do

26              if (mp2 > ndo2) then
                    exit
                end if

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

                        if (mlat == 0) then
                            exit !go to 28
                        end if

                        ve(imid, 2*mp1-2, k) = ve(imid, 2*mp1-2, k) &
                            +br(mp1, np1, k)*vb(imid, mn)
                        ve(imid, 2*mp1-1, k) = ve(imid, 2*mp1-1, k) &
                            +bi(mp1, np1, k)*vb(imid, mn)
                        we(imid, 2*mp1-2, k) = we(imid, 2*mp1-2, k) &
                            -cr(mp1, np1, k)*vb(imid, mn)
                        we(imid, 2*mp1-1, k) = we(imid, 2*mp1-1, k) &
                            -ci(mp1, np1, k)*vb(imid, mn)
                    end do
                end do
            end do
        case(1)
            !
            !==> case ityp=1   no symmetries,  cr and ci equal zero
            !
            !    case m = 0
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
            if (mmax < 2) then
                exit vector_symmetries
            end if

            do mp1=2, mmax
                m = mp1-1
                !     mb = m*(nlat-1)-(m*(m-1))/2
                mb = m*nlat-(m*(m+1))/2
                mp2 = mp1+1

                if (mp1 > ndo1) then
                    go to 126
                end if

                do k=1, nt
                    do np1=mp1, ndo1, 2
                        mn = mb+np1
                        do i=1, imm1
                            vo(i, 2*mp1-2, k) = vo(i, 2*mp1-2, k)+br(mp1, np1, k)*vb(i, mn)
                            vo(i, 2*mp1-1, k) = vo(i, 2*mp1-1, k)+bi(mp1, np1, k)*vb(i, mn)
                            we(i, 2*mp1-2, k) = we(i, 2*mp1-2, k)-bi(mp1, np1, k)*wb(i, mn)
                            we(i, 2*mp1-1, k) = we(i, 2*mp1-1, k)+br(mp1, np1, k)*wb(i, mn)
                        end do

                        if (mlat == 0) then
                            exit
                        end if

                        we(imid, 2*mp1-2, k) = we(imid, 2*mp1-2, k) &
                            -bi(mp1, np1, k)*wb(imid, mn)
                        we(imid, 2*mp1-1, k) = we(imid, 2*mp1-1, k) &
                            +br(mp1, np1, k)*wb(imid, mn)
                    end do
                end do

126             if (mp2 > ndo2) then
                    exit
                end if

                do k=1, nt
                    do np1=mp2, ndo2, 2
                        mn = mb+np1
                        do i=1, imm1
                            ve(i, 2*mp1-2, k) = ve(i, 2*mp1-2, k)+br(mp1, np1, k)*vb(i, mn)
                            ve(i, 2*mp1-1, k) = ve(i, 2*mp1-1, k)+bi(mp1, np1, k)*vb(i, mn)
                            wo(i, 2*mp1-2, k) = wo(i, 2*mp1-2, k)-bi(mp1, np1, k)*wb(i, mn)
                            wo(i, 2*mp1-1, k) = wo(i, 2*mp1-1, k)+br(mp1, np1, k)*wb(i, mn)
                        end do

                        if (mlat == 0) then
                            exit
                        end if

                        ve(imid, 2*mp1-2, k) = ve(imid, 2*mp1-2, k) &
                            +br(mp1, np1, k)*vb(imid, mn)
                        ve(imid, 2*mp1-1, k) = ve(imid, 2*mp1-1, k) &
                            +bi(mp1, np1, k)*vb(imid, mn)
                    end do
                end do
            end do
        case(2)
            !
            !==> case ityp=2   no symmetries,  br and bi are equal to zero
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
            !==> case m = 1 through nlat-1
            !
            if (mmax < 2) then
                exit vector_symmetries
            end if

            do mp1=2, mmax
                m = mp1-1
                !     mb = m*(nlat-1)-(m*(m-1))/2
                mb = m*nlat-(m*(m+1))/2
                mp2 = mp1+1

                if (mp1 > ndo1) then
                    go to 226
                end if

                do k=1, nt
                    do np1=mp1, ndo1, 2
                        mn = mb+np1
                        do i=1, imm1
                            ve(i, 2*mp1-2, k) = ve(i, 2*mp1-2, k)-ci(mp1, np1, k)*wb(i, mn)
                            ve(i, 2*mp1-1, k) = ve(i, 2*mp1-1, k)+cr(mp1, np1, k)*wb(i, mn)
                            wo(i, 2*mp1-2, k) = wo(i, 2*mp1-2, k)-cr(mp1, np1, k)*vb(i, mn)
                            wo(i, 2*mp1-1, k) = wo(i, 2*mp1-1, k)-ci(mp1, np1, k)*vb(i, mn)
                        end do

                        if (mlat == 0) then
                            exit
                        end if

                        ve(imid, 2*mp1-2, k) = ve(imid, 2*mp1-2, k) &
                            -ci(mp1, np1, k)*wb(imid, mn)
                        ve(imid, 2*mp1-1, k) = ve(imid, 2*mp1-1, k) &
                            +cr(mp1, np1, k)*wb(imid, mn)
                    end do
                end do

226             if (mp2 > ndo2) then
                    exit
                end if

                do k=1, nt
                    do np1=mp2, ndo2, 2
                        mn = mb+np1
                        do i=1, imm1
                            vo(i, 2*mp1-2, k) = vo(i, 2*mp1-2, k)-ci(mp1, np1, k)*wb(i, mn)
                            vo(i, 2*mp1-1, k) = vo(i, 2*mp1-1, k)+cr(mp1, np1, k)*wb(i, mn)
                            we(i, 2*mp1-2, k) = we(i, 2*mp1-2, k)-cr(mp1, np1, k)*vb(i, mn)
                            we(i, 2*mp1-1, k) = we(i, 2*mp1-1, k)-ci(mp1, np1, k)*vb(i, mn)
                        end do

                        if (mlat == 0) then
                            exit !go to 228
                        end if

                        we(imid, 2*mp1-2, k) = we(imid, 2*mp1-2, k) &
                            -cr(mp1, np1, k)*vb(imid, mn)
                        we(imid, 2*mp1-1, k) = we(imid, 2*mp1-1, k) &
                            -ci(mp1, np1, k)*vb(imid, mn)
                    end do
                end do
            end do
        case (3)
            !
            !==> case ityp=3   v even,  w odd
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
            !==> case m = 1 through nlat-1
            !
            if (mmax < 2) then
                exit vector_symmetries
            end if

            do mp1=2, mmax
                m = mp1-1
                !     mb = m*(nlat-1)-(m*(m-1))/2
                mb = m*nlat-(m*(m+1))/2
                mp2 = mp1+1

                if (mp1 > ndo1) then
                    go to 326
                end if

                do k=1, nt
                    do np1=mp1, ndo1, 2
                        mn = mb+np1
                        do i=1, imm1
                            ve(i, 2*mp1-2, k) = ve(i, 2*mp1-2, k)-ci(mp1, np1, k)*wb(i, mn)
                            ve(i, 2*mp1-1, k) = ve(i, 2*mp1-1, k)+cr(mp1, np1, k)*wb(i, mn)
                            wo(i, 2*mp1-2, k) = wo(i, 2*mp1-2, k)-cr(mp1, np1, k)*vb(i, mn)
                            wo(i, 2*mp1-1, k) = wo(i, 2*mp1-1, k)-ci(mp1, np1, k)*vb(i, mn)
                        end do

                        if (mlat == 0) then
                            exit
                        end if

                        ve(imid, 2*mp1-2, k) = ve(imid, 2*mp1-2, k) &
                            -ci(mp1, np1, k)*wb(imid, mn)
                        ve(imid, 2*mp1-1, k) = ve(imid, 2*mp1-1, k) &
                            +cr(mp1, np1, k)*wb(imid, mn)
                    end do
                end do

326             if (mp2 > ndo2) then
                    exit
                end if

                do k=1, nt
                    do np1=mp2, ndo2, 2
                        mn = mb+np1
                        do i=1, imm1
                            ve(i, 2*mp1-2, k) = ve(i, 2*mp1-2, k)+br(mp1, np1, k)*vb(i, mn)
                            ve(i, 2*mp1-1, k) = ve(i, 2*mp1-1, k)+bi(mp1, np1, k)*vb(i, mn)
                            wo(i, 2*mp1-2, k) = wo(i, 2*mp1-2, k)-bi(mp1, np1, k)*wb(i, mn)
                            wo(i, 2*mp1-1, k) = wo(i, 2*mp1-1, k)+br(mp1, np1, k)*wb(i, mn)
                        end do

                        if (mlat == 0) then
                            exit !go to 328
                        end if

                        ve(imid, 2*mp1-2, k) = ve(imid, 2*mp1-2, k) &
                            +br(mp1, np1, k)*vb(imid, mn)
                        ve(imid, 2*mp1-1, k) = ve(imid, 2*mp1-1, k) &
                            +bi(mp1, np1, k)*vb(imid, mn)
                    end do
                end do
            end do
        case (4)
            !
            !==> case ityp=4   v even,  w odd, and both cr and ci equal zero
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
            if (mmax < 2) then
                exit vector_symmetries
            end if

            do mp1=2, mmax
                m = mp1-1
                !     mb = m*(nlat-1)-(m*(m-1))/2
                mb = m*nlat-(m*(m+1))/2
                mp2 = mp1+1

                if (mp2 > ndo2) then
                    exit !go to 430
                end if

                do k=1, nt
                    do np1=mp2, ndo2, 2
                        mn = mb+np1
                        do i=1, imm1
                            ve(i, 2*mp1-2, k) = ve(i, 2*mp1-2, k)+br(mp1, np1, k)*vb(i, mn)
                            ve(i, 2*mp1-1, k) = ve(i, 2*mp1-1, k)+bi(mp1, np1, k)*vb(i, mn)
                            wo(i, 2*mp1-2, k) = wo(i, 2*mp1-2, k)-bi(mp1, np1, k)*wb(i, mn)
                            wo(i, 2*mp1-1, k) = wo(i, 2*mp1-1, k)+br(mp1, np1, k)*wb(i, mn)
                        end do

                        if (mlat == 0) then
                            exit
                        end if

                        ve(imid, 2*mp1-2, k) = ve(imid, 2*mp1-2, k) &
                            +br(mp1, np1, k)*vb(imid, mn)
                        ve(imid, 2*mp1-1, k) = ve(imid, 2*mp1-1, k) &
                            +bi(mp1, np1, k)*vb(imid, mn)
                    end do
                end do
            end do
        case (5)
            !
            !==> case ityp=5   v even,  w odd,     br and bi equal zero
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
            if (mmax < 2) then
                exit vector_symmetries
            end if

            do mp1=2, mmax
                m = mp1-1
                !     mb = m*(nlat-1)-(m*(m-1))/2
                mb = m*nlat-(m*(m+1))/2
                mp2 = mp1+1

                if (mp1 > ndo1) then
                    exit
                end if

                do k=1, nt
                    do np1=mp1, ndo1, 2
                        mn = mb+np1
                        do i=1, imm1
                            ve(i, 2*mp1-2, k) = ve(i, 2*mp1-2, k)-ci(mp1, np1, k)*wb(i, mn)
                            ve(i, 2*mp1-1, k) = ve(i, 2*mp1-1, k)+cr(mp1, np1, k)*wb(i, mn)
                            wo(i, 2*mp1-2, k) = wo(i, 2*mp1-2, k)-cr(mp1, np1, k)*vb(i, mn)
                            wo(i, 2*mp1-1, k) = wo(i, 2*mp1-1, k)-ci(mp1, np1, k)*vb(i, mn)
                        end do

                        if (mlat == 0) then
                            exit
                        end if

                        ve(imid, 2*mp1-2, k) = ve(imid, 2*mp1-2, k) &
                            -ci(mp1, np1, k)*wb(imid, mn)
                        ve(imid, 2*mp1-1, k) = ve(imid, 2*mp1-1, k) &
                            +cr(mp1, np1, k)*wb(imid, mn)
                    end do
                end do
            end do
        case (6)
            !
            !==> case ityp=6   v odd  ,  w even
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
            if (mmax < 2) then
                exit vector_symmetries
            end if

            do mp1=2, mmax
                m = mp1-1
                !     mb = m*(nlat-1)-(m*(m-1))/2
                mb = m*nlat-(m*(m+1))/2
                mp2 = mp1+1

                if (mp1 > ndo1) then
                    go to 626
                end if

                do k=1, nt
                    do np1=mp1, ndo1, 2
                        mn = mb+np1
                        do i=1, imm1
                            vo(i, 2*mp1-2, k) = vo(i, 2*mp1-2, k)+br(mp1, np1, k)*vb(i, mn)
                            vo(i, 2*mp1-1, k) = vo(i, 2*mp1-1, k)+bi(mp1, np1, k)*vb(i, mn)
                            we(i, 2*mp1-2, k) = we(i, 2*mp1-2, k)-bi(mp1, np1, k)*wb(i, mn)
                            we(i, 2*mp1-1, k) = we(i, 2*mp1-1, k)+br(mp1, np1, k)*wb(i, mn)
                        end do

                        if (mlat == 0) then
                            exit!go to 624
                        end if

                        we(imid, 2*mp1-2, k) = we(imid, 2*mp1-2, k) &
                            -bi(mp1, np1, k)*wb(imid, mn)
                        we(imid, 2*mp1-1, k) = we(imid, 2*mp1-1, k) &
                            +br(mp1, np1, k)*wb(imid, mn)
                    end do
                end do

626             if (mp2 > ndo2) then
                    exit
                end if


                do k=1, nt
                    do np1=mp2, ndo2, 2
                        mn = mb+np1
                        do i=1, imm1
                            vo(i, 2*mp1-2, k) = vo(i, 2*mp1-2, k)-ci(mp1, np1, k)*wb(i, mn)
                            vo(i, 2*mp1-1, k) = vo(i, 2*mp1-1, k)+cr(mp1, np1, k)*wb(i, mn)
                            we(i, 2*mp1-2, k) = we(i, 2*mp1-2, k)-cr(mp1, np1, k)*vb(i, mn)
                            we(i, 2*mp1-1, k) = we(i, 2*mp1-1, k)-ci(mp1, np1, k)*vb(i, mn)
                        end do

                        if (mlat == 0) then
                            exit
                        end if

                        we(imid, 2*mp1-2, k) = we(imid, 2*mp1-2, k) &
                            -cr(mp1, np1, k)*vb(imid, mn)
                        we(imid, 2*mp1-1, k) = we(imid, 2*mp1-1, k) &
                            -ci(mp1, np1, k)*vb(imid, mn)
                    end do
                end do
            end do
        case(7)
            !
            !==> case ityp=7   v odd, w even   cr and ci equal zero
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
            !==> case m = 1 through nlat-1
            !
            if (mmax < 2) then
                exit vector_symmetries
            end if

            do mp1=2, mmax
                m = mp1-1
                !     mb = m*(nlat-1)-(m*(m-1))/2
                mb = m*nlat-(m*(m+1))/2
                mp2 = mp1+1

                if (mp1 > ndo1) then
                    exit
                end if

                do k=1, nt
                    do np1=mp1, ndo1, 2
                        mn = mb+np1
                        do i=1, imm1
                            vo(i, 2*mp1-2, k) = vo(i, 2*mp1-2, k)+br(mp1, np1, k)*vb(i, mn)
                            vo(i, 2*mp1-1, k) = vo(i, 2*mp1-1, k)+bi(mp1, np1, k)*vb(i, mn)
                            we(i, 2*mp1-2, k) = we(i, 2*mp1-2, k)-bi(mp1, np1, k)*wb(i, mn)
                            we(i, 2*mp1-1, k) = we(i, 2*mp1-1, k)+br(mp1, np1, k)*wb(i, mn)
                        end do

                        if (mlat == 0) then
                            exit
                        end if

                        we(imid, 2*mp1-2, k) = we(imid, 2*mp1-2, k) &
                            -bi(mp1, np1, k)*wb(imid, mn)
                        we(imid, 2*mp1-1, k) = we(imid, 2*mp1-1, k) &
                            +br(mp1, np1, k)*wb(imid, mn)
                    end do
                end do
            end do
        case (8)
            !
            !==> case ityp=8   v odd,  w even   br and bi equal zero
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
            if (mmax < 2) then
                exit vector_symmetries
            end if

            do mp1=2, mmax
                m = mp1-1
                !     mb = m*(nlat-1)-(m*(m-1))/2
                mb = m*nlat-(m*(m+1))/2
                mp2 = mp1+1

                if (mp2 > ndo2) then
                    exit
                end if

                do k=1, nt
                    do np1=mp2, ndo2, 2
                        mn = mb+np1
                        do i=1, imm1
                            vo(i, 2*mp1-2, k) = vo(i, 2*mp1-2, k)-ci(mp1, np1, k)*wb(i, mn)
                            vo(i, 2*mp1-1, k) = vo(i, 2*mp1-1, k)+cr(mp1, np1, k)*wb(i, mn)
                            we(i, 2*mp1-2, k) = we(i, 2*mp1-2, k)-cr(mp1, np1, k)*vb(i, mn)
                            we(i, 2*mp1-1, k) = we(i, 2*mp1-1, k)-ci(mp1, np1, k)*vb(i, mn)
                        end do

                        if (mlat == 0) then
                            exit
                        end if

                        we(imid, 2*mp1-2, k) = we(imid, 2*mp1-2, k) &
                            -cr(mp1, np1, k)*vb(imid, mn)
                        we(imid, 2*mp1-1, k) = we(imid, 2*mp1-1, k) &
                            -ci(mp1, np1, k)*vb(imid, mn)
                    end do
                end do
            end do
    end select vector_symmetries

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
        go to 13
    end if

    do k=1, nt
        do j=1, nlon
            do i=1, imm1
                v(i, j, k) = .5*ve(i, j, k)
                w(i, j, k) = .5*we(i, j, k)
            end do
        end do
    end do

13  if (mlat == 0) then
        return
    end if

    do k=1, nt
        do j=1, nlon
            v(imid, j, k) = .5*ve(imid, j, k)
            w(imid, j, k) = .5*we(imid, j, k)
        end do
    end do

end subroutine vhsgs1



subroutine vhsgsi(nlat, nlon, wvhsgs, lvhsgs, dwork, ldwork, ierror)
    !
    !     subroutine vhsfsi computes the gaussian points theta, gauss
    !     weights wts, and the components vb and wb of the vector
    !     harmonics. all quantities are computed internally in double
    !     precision but returned in single precision and are therfore
    !     accurate to single precision.
    !
    !     set imid = (nlat+1)/2 and lmn=(nlat*(nlat+1))/2 then
    !     wvhsgs must have 2*(imid*lmn+nlat)+nlon+15 locations
    !
    !     real array dwork must have
    !       3*nlat*(nlat+1)+5*nlat+1 = nlat*(3*nlat+8)+1
    !     locations which is determined by the size of dthet,
    !     dwts, dwork, and dpbar in vhsgs1
    !
    dimension wvhsgs(*)
    real dwork(*)
    ierror = 1
    if (nlat < 3) return
    ierror = 2
    if (nlon < 1) return
    ierror = 3
    imid = (nlat+1)/2
    lmn = (nlat*(nlat+1))/2
    if (lvhsgs < 2*(imid*lmn)+nlon+15) return
    ierror = 4
    if (ldwork < (nlat*3*(nlat+3)+2)/2) return
    ierror = 0
    !
    !     set saved work space pointers
    !
    jw1 = 1
    jw2 = jw1+imid*lmn
    jw3 = jw2+imid*lmn
    !
    !     set unsaved work space pointers
    !
    iw1 = 1
    iw2 = iw1+nlat
    iw3 = iw2+nlat
    iw4 = iw3+3*imid*nlat
    !     iw2 = iw1+nlat+nlat
    !     iw3 = iw2+nlat+nlat
    !     iw4 = iw3+6*imid*nlat
    call vhgsi1(nlat, imid, wvhsgs(jw1), wvhsgs(jw2), &
        dwork(iw1), dwork(iw2), dwork(iw3), dwork(iw4))
    call hrffti(nlon, wvhsgs(jw3))

end subroutine vhsgsi



subroutine vhgsi1(nlat, imid, vb, wb, dthet, dwts, dpbar, work)
    dimension vb(imid, *), wb(imid, *)
    real abel, bbel, cbel, ssqr2, dcf
    real dthet(*), dwts(*), dpbar(imid, nlat, 3), work(*)
    !
    !     compute gauss points and weights
    !     use dpbar (length 3*nnlat*(nnlat+1)) as work space for gaqd
    !
    lwk = nlat*(nlat+2)
    call gaqd(nlat, dthet, dwts, dpbar, lwk, ierror)
    !
    !     compute associated legendre functions
    !
    !     compute m=n=0 legendre polynomials for all theta(i)
    !
    ssqr2 = 1.0/sqrt(2.0)
    dpbar(1:imid, 1, 1) = ssqr2
    vb(1:imid, 1) = 0.0
    wb(1:imid, 1) = 0.0
    !
    !     main loop for remaining vb, and wb
    !
    do n=1, nlat-1
        nm = mod(n-2, 3)+1
        nz = mod(n-1, 3)+1
        np = mod(n, 3)+1
        !
        !==> compute dpbar for m=0
        !
        call dnlfk(0, n, work)
        mn = indx(0, n, nlat)
        do i=1, imid
            call dnlft(0, n, dthet(i), work, dpbar(i, 1, np))
        end do
        !
        !==> compute dpbar for m=1
        !
        call dnlfk(1, n, work)
        mn = indx(1, n, nlat)
        do i=1, imid
            call dnlft(1, n, dthet(i), work, dpbar(i, 2, np))
        end do
        !
        !==> compute and store dpbar for m=2, n
        !
        if (n >= 2) then
            do m=2, n
                abel = sqrt(real((2*n+1)*(m+n-2)*(m+n-3))/ &
                    real((2*n-3)*(m+n-1)*(m+n)))
                bbel = sqrt(real((2*n+1)*(n-m-1)*(n-m))/ &
                    real((2*n-3)*(m+n-1)*(m+n)))
                cbel = sqrt(real((n-m+1)*(n-m+2))/ &
                    real((m+n-1)*(m+n)))
                id = indx(m, n, nlat)
                if (m >= n-1) then
                    dpbar(1:imid, m+1, np) = &
                        abel*dpbar(1:imid, m-1, nm)-cbel*dpbar(1:imid, m-1, np)
                else
                    dpbar(1:imid, m+1, np) = &
                        abel*dpbar(1:imid, m-1, nm)+bbel*dpbar(1:imid, m+1, nm) &
                        -cbel*dpbar(1:imid, m-1, np)
                end if
            end do
        end if
        !
        !==> compute the derivative of the functions
        !
        ix = indx(0, n, nlat)
        iy = indx(n, n, nlat)
        vb(1:imid, ix) = -dpbar(1:imid, 2, np)
        vb(1:imid, iy) = dpbar(1:imid, n, np)/sqrt(real(2*(n+1)))

        if (n /= 1) then
            dcf = sqrt(real(4*n*(n+1)))
            do m=1, n-1
                ix = indx(m, n, nlat)
                abel = sqrt(real((n+m)*(n-m+1)))/dcf
                bbel = sqrt(real((n-m)*(n+m+1)))/dcf
                vb(1:imid, ix) = &
                    abel*dpbar(1:imid, m, np)-bbel*dpbar(1:imid, m+2, np)
            end do
        end if
        !
        !==> compute the vector harmonic w(theta) = m*pbar/cos(theta)
        !
        !     set wb=0 for m=0
        !
        ix = indx(0, n, nlat)
        wb(1:imid, ix) = 0.0
        !
        !==> compute wb for m=1, n
        !
        dcf = sqrt(real(n+n+1)/real(4*n*(n+1)*(n+n-1)))
        do m=1, n
            ix = indx(m, n, nlat)
            abel = dcf*sqrt(real((n+m)*(n+m-1)))
            bbel = dcf*sqrt(real((n-m)*(n-m-1)))
            if (m >= n-1) then
                wb(1:imid, ix) = abel*dpbar(1:imid, m, nz)
            else
                wb(1:imid, ix) = &
                    abel*dpbar(1:imid, m, nz) + bbel*dpbar(1:imid, m+2, nz)
            end if
        end do
    end do

end subroutine vhgsi1
