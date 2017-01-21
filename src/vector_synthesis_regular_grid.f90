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
! ... file vhsec.f
!
!     this file contains code and documentation for subroutines
!     vhsec and vhseci
!
! ... files which must be loaded with vhsec.f
!
!     type_SpherepackAux.f, type_HFFTpack.f
!   
!     subroutine vhsec(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, 
!    +                 mdab, ndab, wvhsec, lvhsec, work, lwork, ierror)
!
!     subroutine vhsec performs the vector spherical harmonic synthesis
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
!            nlat/2 if nlat is even or(nlat+1)/2 if nlat is odd.
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
!                 i=1, ...,(nlat+1)/2 and j=1, ..., nlon. if nlat is
!                 even the synthesis is performed on the the arrays
!                 v(i, j), w(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
!
!            = 4  v is symmetric and w is antisymmetric about the 
!                 equator. the synthesis is performed on the northern
!                 hemisphere only.  i.e., if nlat is odd the synthesis
!                 is performed on the arrays v(i, j), w(i, j) for 
!                 i=1, ...,(nlat+1)/2 and j=1, ..., nlon. if nlat is
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
!                 i=1, ...,(nlat+1)/2 and j=1, ..., nlon. if nlat is
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
!                 i=1, ...,(nlat+1)/2 and j=1, ..., nlon. if nlat is
!                 even the synthesis is performed on the the arrays
!                 v(i, j), w(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
!
!            = 7  v is antisymmetric and w is symmetric about the 
!                 equator. the synthesis is performed on the northern
!                 hemisphere only.  i.e., if nlat is odd the synthesis
!                 is performed on the arrays v(i, j), w(i, j) for 
!                 i=1, ...,(nlat+1)/2 and j=1, ..., nlon. if nlat is
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
!                 i=1, ...,(nlat+1)/2 and j=1, ..., nlon. if nlat is
!                 even the synthesis is performed on the the arrays
!                 v(i, j), w(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
!                 the divergence of (v, w) is zero. i.e., 
!                 (d/dtheta (sin(theta) v) + dw/dphi)/sin(theta) = 0. 
!                 the coefficients br and bi are zero.
!
!
!     nt     the number of syntheses.  in the program that calls vhsec, 
!            the arrays v, w, br, bi, cr, and ci can be three dimensional
!            in which case multiple syntheses will be performed.
!            the third index is the synthesis index which assumes the 
!            values k=1, ..., nt.  for a single synthesis set nt=1. the
!            discription of the remaining parameters is simplified
!            by assuming that nt=1 or that all the arrays are two
!            dimensional.
!
!     idvw   the first dimension of the arrays v, w as it appears in
!            the program that calls vhsec. if ityp .le. 2 then idvw
!            must be at least nlat.  if ityp .gt. 2 and nlat is
!            even then idvw must be at least nlat/2. if ityp .gt. 2
!            and nlat is odd then idvw must be at least(nlat+1)/2.
!
!     jdvw   the second dimension of the arrays v, w as it appears in
!            the program that calls vhsec. jdvw must be at least nlon.
!
!     br, bi  two or three dimensional arrays (see input parameter nt)
!     cr, ci  that contain the vector spherical harmonic coefficients
!            in the spectral representation of v(i, j) and w(i, j) given
!            below at the discription of output parameters v and w.
!
!     mdab   the first dimension of the arrays br, bi, cr, and ci as it
!            appears in the program that calls vhsec. mdab must be at
!            least min(nlat, nlon/2) if nlon is even or at least
!            min(nlat, (nlon+1)/2) if nlon is odd.
!
!     ndab   the second dimension of the arrays br, bi, cr, and ci as it
!            appears in the program that calls vhsec. ndab must be at
!            least nlat.
!
!     wvhsec an array which must be initialized by subroutine vhseci.
!            once initialized, wvhsec can be used repeatedly by vhsec
!            as long as nlon and nlat remain unchanged.  wvhsec must
!            not be altered between calls of vhsec.
!
!     lvhsec the dimension of the array wvhsec as it appears in the
!            program that calls vhsec. define
!
!               l1 = min(nlat, nlon/2) if nlon is even or
!               l1 = min(nlat, (nlon+1)/2) if nlon is odd
!
!            and
!
!               l2 = nlat/2        if nlat is even or
!               l2 =(nlat+1)/2    if nlat is odd
!
!            then lvhsec must be at least
!
!               4*nlat*l2+3*max(l1-2, 0)*(nlat+nlat-l1-1)+nlon+15
!
!
!     work   a work array that does not have to be saved.
!
!     lwork  the dimension of the array work as it appears in the
!            program that calls vhsec. define
!
!               l2 = nlat/2        if nlat is even or
!               l2 =(nlat+1)/2    if nlat is odd
!
!            if ityp .le. 2 then lwork must be at least
!
!                    nlat*(2*nt*nlon+max(6*l2, nlon))
!
!            if ityp .gt. 2 then lwork must be at least
!
!                    l2*(2*nt*nlon+max(6*nlat, nlon))
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
!            = 9  error in the specification of lvhsec
!            = 10 error in the specification of lwork
!
!
! *******************************************************************
!
!     subroutine vhseci(nlat, nlon, wvhsec, lvhsec, dwork, ldwork, ierror)
!
!     subroutine vhseci initializes the array wvhsec which can then be
!     used repeatedly by subroutine vhsec until nlat or nlon is changed.
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
!            nlat/2 if nlat is even or(nlat+1)/2 if nlat is odd.
!
!     nlon   the number of distinct londitude points.  nlon determines
!            the grid increment in longitude as 2*pi/nlon. for example
!            nlon = 72 for a five degree grid. nlon must be greater
!            than zero. the axisymmetric case corresponds to nlon=1.
!            the efficiency of the computation is improved when nlon
!            is a product of small prime numbers.
!
!     lvhsec the dimension of the array wvhsec as it appears in the
!            program that calls vhsec. define
!
!               l1 = min(nlat, nlon/2) if nlon is even or
!               l1 = min(nlat, (nlon+1)/2) if nlon is odd
!
!            and
!
!               l2 = nlat/2        if nlat is even or
!               l2 =(nlat+1)/2    if nlat is odd
!
!            then lvhsec must be at least
!
!            4*nlat*l2+3*max(l1-2, 0)*(nlat+nlat-l1-1)+nlon+15
!
!
!     dwork  a real work array that does not have to be saved.
!
!     ldwork the dimension of the array dwork as it appears in the
!            program that calls vhsec. ldwork must be at least
!            2*(nlat+2)
!
!     **************************************************************
!
!     output parameters
!
!     wvhsec an array which is initialized for use by subroutine vhsec.
!            once initialized, wvhsec can be used repeatedly by vhsec
!            as long as nlat or nlon remain unchanged.  wvhsec must not
!            be altered between calls of vhsec.
!
!
!     ierror = 0  no errors
!            = 1  error in the specification of nlat
!            = 2  error in the specification of nlon
!            = 3  error in the specification of lvhsec
!            = 4  error in the specification of ldwork
!
!
!
submodule(vector_synthesis_routines) vector_synthesis_regular_grid

contains

    module subroutine vhsec(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, &
        mdab, ndab, wvhsec, lvhsec, work, lwork, ierror)

        ! Dummy arguments
        real(wp) :: br(mdab, ndab, nt), bi(mdab, ndab, nt)
        real(wp) :: cr(mdab, ndab, nt), ci(mdab, ndab, nt)
        integer(ip) :: idvw
        integer(ip) :: ierror
        integer(ip) :: ityp
        integer(ip) :: lvhsec
        integer(ip) :: lwork
        integer(ip) :: mdab
        integer(ip) :: ndab
        integer(ip) :: nlat
        integer(ip) :: nlon
        integer(ip) :: nt
        real(wp) :: v(idvw, jdvw, nt), w(idvw, jdvw, nt)
        real(wp) :: work(lwork), wvhsec(lvhsec)

        ! Local variables
        integer(ip) :: idv
        integer(ip) :: imid
        integer(ip) :: ist
        integer(ip) :: iw1
        integer(ip) :: iw2
        integer(ip) :: iw3
        integer(ip) :: iw4
        integer(ip) :: iw5
        integer(ip) :: jdvw
        integer(ip) :: jw1
        integer(ip) :: jw2
        integer(ip) :: labc
        integer(ip) :: lnl
        integer(ip) :: lwzvin
        integer(ip) :: lzz1
        integer(ip) :: mmax

        ierror = 1
        if(nlat < 3) return
        ierror = 2
        if (nlon < 1) return
        ierror = 3
        if (ityp<0 .or. ityp>8) return
        ierror = 4
        if (nt < 0) return
        ierror = 5
        imid =(nlat+1)/2
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
        lzz1 = 2*nlat*imid
        labc = 3*(max(mmax-2, 0)*(nlat+nlat-mmax-1))/2
        if (lvhsec < 2*(lzz1+labc)+nlon+15) return
        ierror = 10
        if (ityp <= 2 .and. &
            lwork < nlat*(2*nt*nlon+max(6*imid, nlon))) return
        if (ityp > 2 .and. &
            lwork < imid*(2*nt*nlon+max(6*nlat, nlon))) return
        ierror = 0
        idv = nlat
        if (ityp > 2) idv = imid
        lnl = nt*idv*nlon
        ist = 0
        if (ityp <= 2) ist = imid
        iw1 = ist+1
        iw2 = lnl+1
        iw3 = iw2+ist
        iw4 = iw2+lnl
        iw5 = iw4+3*imid*nlat
        lzz1 = 2*nlat*imid
        labc = 3*(max(mmax-2, 0)*(2*nlat-mmax-1))/2
        lwzvin = lzz1+labc
        jw1 = lwzvin+1
        jw2 = jw1+lwzvin

        call vhsec_lower_routine(nlat, nlon, ityp, nt, imid, idvw, jdvw, v, w, mdab, ndab, &
            br, bi, cr, ci, idv, work, work(iw1), work(iw2), work(iw3), &
            work(iw4), work(iw5), wvhsec, wvhsec(jw1), wvhsec(jw2))

    end subroutine vhsec

    module subroutine vhseci(nlat, nlon, wvhsec, lvhsec, dwork, ldwork, ierror)

        ! Dummy arguments
        integer(ip) :: ierror
        integer(ip) :: nlat
        integer(ip) :: nlon
        real(wp) :: wvhsec(lvhsec)
        real(wp) :: dwork(ldwork)
        integer(ip) :: ldwork
        integer(ip) :: lvhsec

        ! Local variables
        integer(ip) :: imid
        integer(ip) :: iw1
        integer(ip) :: iw2
        integer(ip) :: labc
        integer(ip) :: lwvbin
        integer(ip) :: lzz1
        integer(ip) :: mmax
        type(SpherepackAux) :: sphere_aux

        imid =(nlat+1)/2
        lzz1 = 2*nlat*imid
        mmax = min(nlat, (nlon+1)/2)
        labc = 3*(max(mmax-2, 0)*(2*nlat-mmax-1))/2

        ierror = 1
        if(nlat < 3) return
        ierror = 2
        if (nlon < 1) return
        ierror = 3
        if (lvhsec < 2*(lzz1+labc)+nlon+15) return
        ierror = 4
        if (ldwork < 2*nlat+2) return
        ierror = 0

        !  Set workspace pointers
        lwvbin = lzz1+labc
        iw1 = lwvbin+1
        iw2 = iw1+lwvbin

        call sphere_aux%vbinit(nlat, nlon, wvhsec, dwork)

        call sphere_aux%wbinit(nlat, nlon, wvhsec(iw1), dwork)

        call sphere_aux%hfft%initialize(nlon, wvhsec(iw2))

    end subroutine vhseci

    subroutine vhsec_lower_routine(nlat, nlon, ityp, nt, imid, idvw, jdvw, v, w, mdab, &
        ndab, br, bi, cr, ci, idv, ve, vo, we, wo, vb, wb, wvbin, wwbin, wrfft)

        real(wp) :: bi
        real(wp) :: br
        real(wp) :: ci
        real(wp) :: cr
        integer(ip) :: i
        integer(ip) :: idv
        integer(ip) :: idvw
        integer(ip) :: imid
        integer(ip) :: imm1
        integer(ip) :: ityp
        integer(ip) :: itypp
        integer(ip) :: iv
        integer(ip) :: iw
        integer(ip) :: j
        integer(ip) :: jdvw
        integer(ip) :: k
        integer(ip) :: m
        integer(ip) :: mdab
        integer(ip) :: mlat
        integer(ip) :: mlon
        integer(ip) :: mmax
        integer(ip) :: mp1
        integer(ip) :: mp2
        integer(ip) :: ndab
        integer(ip) :: ndo1
        integer(ip) :: ndo2
        integer(ip) :: nlat
        integer(ip) :: nlon
        integer(ip) :: nlp1
        integer(ip) :: np1
        integer(ip) :: nt
        real(wp) :: v
        real(wp) :: vb
        real(wp) :: ve
        real(wp) :: vo
        real(wp) :: w
        real(wp) :: wb
        real(wp) :: we
        real(wp) :: wo
        real(wp) :: wrfft
        real(wp) :: wvbin
        real(wp) :: wwbin
        dimension v(idvw, jdvw, *), w(idvw, jdvw, *), br(mdab, ndab, *), &
            bi(mdab, ndab, *), cr(mdab, ndab, *), ci(mdab, ndab, *), &
            ve(idv, nlon, *), vo(idv, nlon, *), we(idv, nlon, *), &
            wo(idv, nlon, *), wvbin(*), wwbin(*), wrfft(*), &
            vb(imid, nlat, 3), wb(imid, nlat, 3)

        type(SpherepackAux) :: sphere_aux

        nlp1 = nlat+1
        mlat = mod(nlat, 2)
        mlon = mod(nlon, 2)
        mmax = min(nlat, (nlon+1)/2)
        imm1 = imid
        if (mlat /= 0) imm1 = imid-1
        do k=1, nt
            do j=1, nlon
                do i=1, idv
                    ve(i, j, k) = ZERO
                    we(i, j, k) = ZERO
10              continue
                end do
            end do
        end do
        ndo1 = nlat
        ndo2 = nlat
        if (mlat /= 0) ndo1 = nlat-1
        if (mlat == 0) ndo2 = nlat-1
18      itypp = ityp+1
        goto (1, 100, 200, 300, 400, 500, 600, 700, 800), itypp
        !
        !     case ityp=0   no symmetries
        !
1       call sphere_aux%vbin(0, nlat, nlon, 0, vb, iv, wvbin)
        !
        !     case m = 0
        !
        do k=1, nt
            do np1=2, ndo2, 2
                do i=1, imid
                    ve(i, 1, k)=ve(i, 1, k)+br(1, np1, k)*vb(i, np1, iv)
                    we(i, 1, k)=we(i, 1, k)-cr(1, np1, k)*vb(i, np1, iv)
15              continue
                end do
            end do
        end do
        do k=1, nt
            do np1=3, ndo1, 2
                do i=1, imm1
                    vo(i, 1, k)=vo(i, 1, k)+br(1, np1, k)*vb(i, np1, iv)
                    wo(i, 1, k)=wo(i, 1, k)-cr(1, np1, k)*vb(i, np1, iv)
16              continue
                end do
            end do
        end do
        !
        !     case m = 1 through nlat-1
        !
        if (mmax < 2) goto 950
        do mp1=2, mmax
            m = mp1-1
            mp2 = mp1+1
            call sphere_aux%vbin(0, nlat, nlon, m, vb, iv, wvbin)
            call sphere_aux%wbin(0, nlat, nlon, m, wb, iw, wwbin)
            if (mp1 > ndo1) goto 26
            do k=1, nt
                do np1=mp1, ndo1, 2
                    do i=1, imm1
                        vo(i, 2*mp1-2, k) = vo(i, 2*mp1-2, k)+br(mp1, np1, k)*vb(i, np1, iv)
                        ve(i, 2*mp1-2, k) = ve(i, 2*mp1-2, k)-ci(mp1, np1, k)*wb(i, np1, iw)
                        vo(i, 2*mp1-1, k) = vo(i, 2*mp1-1, k)+bi(mp1, np1, k)*vb(i, np1, iv)
                        ve(i, 2*mp1-1, k) = ve(i, 2*mp1-1, k)+cr(mp1, np1, k)*wb(i, np1, iw)
                        wo(i, 2*mp1-2, k) = wo(i, 2*mp1-2, k)-cr(mp1, np1, k)*vb(i, np1, iv)
                        we(i, 2*mp1-2, k) = we(i, 2*mp1-2, k)-bi(mp1, np1, k)*wb(i, np1, iw)
                        wo(i, 2*mp1-1, k) = wo(i, 2*mp1-1, k)-ci(mp1, np1, k)*vb(i, np1, iv)
                        we(i, 2*mp1-1, k) = we(i, 2*mp1-1, k)+br(mp1, np1, k)*wb(i, np1, iw)
23                  continue
                    end do
                    if (mlat == 0) goto 24
                    ve(imid, 2*mp1-2, k) = ve(imid, 2*mp1-2, k) &
                        -ci(mp1, np1, k)*wb(imid, np1, iw)
                    ve(imid, 2*mp1-1, k) = ve(imid, 2*mp1-1, k) &
                        +cr(mp1, np1, k)*wb(imid, np1, iw)
                    we(imid, 2*mp1-2, k) = we(imid, 2*mp1-2, k) &
                        -bi(mp1, np1, k)*wb(imid, np1, iw)
                    we(imid, 2*mp1-1, k) = we(imid, 2*mp1-1, k) &
                        +br(mp1, np1, k)*wb(imid, np1, iw)
24              continue
                end do
25          continue
            end do
26          if (mp2 > ndo2) goto 30
            do k=1, nt
                do np1=mp2, ndo2, 2
                    do i=1, imm1
                        ve(i, 2*mp1-2, k) = ve(i, 2*mp1-2, k)+br(mp1, np1, k)*vb(i, np1, iv)
                        vo(i, 2*mp1-2, k) = vo(i, 2*mp1-2, k)-ci(mp1, np1, k)*wb(i, np1, iw)
                        ve(i, 2*mp1-1, k) = ve(i, 2*mp1-1, k)+bi(mp1, np1, k)*vb(i, np1, iv)
                        vo(i, 2*mp1-1, k) = vo(i, 2*mp1-1, k)+cr(mp1, np1, k)*wb(i, np1, iw)
                        we(i, 2*mp1-2, k) = we(i, 2*mp1-2, k)-cr(mp1, np1, k)*vb(i, np1, iv)
                        wo(i, 2*mp1-2, k) = wo(i, 2*mp1-2, k)-bi(mp1, np1, k)*wb(i, np1, iw)
                        we(i, 2*mp1-1, k) = we(i, 2*mp1-1, k)-ci(mp1, np1, k)*vb(i, np1, iv)
                        wo(i, 2*mp1-1, k) = wo(i, 2*mp1-1, k)+br(mp1, np1, k)*wb(i, np1, iw)
27                  continue
                    end do
                    if (mlat == 0) goto 28
                    ve(imid, 2*mp1-2, k) = ve(imid, 2*mp1-2, k) &
                        +br(mp1, np1, k)*vb(imid, np1, iv)
                    ve(imid, 2*mp1-1, k) = ve(imid, 2*mp1-1, k) &
                        +bi(mp1, np1, k)*vb(imid, np1, iv)
                    we(imid, 2*mp1-2, k) = we(imid, 2*mp1-2, k) &
                        -cr(mp1, np1, k)*vb(imid, np1, iv)
                    we(imid, 2*mp1-1, k) = we(imid, 2*mp1-1, k) &
                        -ci(mp1, np1, k)*vb(imid, np1, iv)
28              continue
                end do
29          continue
            end do
30      continue
        end do
        goto 950
        !
        !     case ityp=1   no symmetries,  cr and ci equal zero
        !
100     call sphere_aux%vbin(0, nlat, nlon, 0, vb, iv, wvbin)
        !
        !     case m = 0
        !
        do k=1, nt
            do np1=2, ndo2, 2
                do i=1, imid
                    ve(i, 1, k)=ve(i, 1, k)+br(1, np1, k)*vb(i, np1, iv)
115             continue
                end do
            end do
        end do
        do k=1, nt
            do np1=3, ndo1, 2
                do i=1, imm1
                    vo(i, 1, k)=vo(i, 1, k)+br(1, np1, k)*vb(i, np1, iv)
116             continue
                end do
            end do
        end do
        !
        !     case m = 1 through nlat-1
        !
        if (mmax < 2) goto 950
        do mp1=2, mmax
            m = mp1-1
            mp2 = mp1+1
            call sphere_aux%vbin(0, nlat, nlon, m, vb, iv, wvbin)
            call sphere_aux%wbin(0, nlat, nlon, m, wb, iw, wwbin)
            if (mp1 > ndo1) goto 126
            do k=1, nt
                do np1=mp1, ndo1, 2
                    do i=1, imm1
                        vo(i, 2*mp1-2, k) = vo(i, 2*mp1-2, k)+br(mp1, np1, k)*vb(i, np1, iv)
                        vo(i, 2*mp1-1, k) = vo(i, 2*mp1-1, k)+bi(mp1, np1, k)*vb(i, np1, iv)
                        we(i, 2*mp1-2, k) = we(i, 2*mp1-2, k)-bi(mp1, np1, k)*wb(i, np1, iw)
                        we(i, 2*mp1-1, k) = we(i, 2*mp1-1, k)+br(mp1, np1, k)*wb(i, np1, iw)
123                 continue
                    end do
                    if (mlat == 0) goto 124
                    we(imid, 2*mp1-2, k) = we(imid, 2*mp1-2, k) &
                        -bi(mp1, np1, k)*wb(imid, np1, iw)
                    we(imid, 2*mp1-1, k) = we(imid, 2*mp1-1, k) &
                        +br(mp1, np1, k)*wb(imid, np1, iw)
124             continue
                end do
125         continue
            end do
126         if (mp2 > ndo2) goto 130
            do k=1, nt
                do np1=mp2, ndo2, 2
                    do i=1, imm1
                        ve(i, 2*mp1-2, k) = ve(i, 2*mp1-2, k)+br(mp1, np1, k)*vb(i, np1, iv)
                        ve(i, 2*mp1-1, k) = ve(i, 2*mp1-1, k)+bi(mp1, np1, k)*vb(i, np1, iv)
                        wo(i, 2*mp1-2, k) = wo(i, 2*mp1-2, k)-bi(mp1, np1, k)*wb(i, np1, iw)
                        wo(i, 2*mp1-1, k) = wo(i, 2*mp1-1, k)+br(mp1, np1, k)*wb(i, np1, iw)
127                 continue
                    end do
                    if (mlat == 0) goto 128
                    ve(imid, 2*mp1-2, k) = ve(imid, 2*mp1-2, k) &
                        +br(mp1, np1, k)*vb(imid, np1, iv)
                    ve(imid, 2*mp1-1, k) = ve(imid, 2*mp1-1, k) &
                        +bi(mp1, np1, k)*vb(imid, np1, iv)
128             continue
                end do
129         continue
            end do
130     continue
        end do
        goto 950
        !
        !     case ityp=2   no symmetries,  br and bi are equal to zero
        !
200     call sphere_aux%vbin(0, nlat, nlon, 0, vb, iv, wvbin)
        !
        !     case m = 0
        !
        do k=1, nt
            do np1=2, ndo2, 2
                do i=1, imid
                    we(i, 1, k)=we(i, 1, k)-cr(1, np1, k)*vb(i, np1, iv)
215             continue
                end do
            end do
        end do
        do k=1, nt
            do np1=3, ndo1, 2
                do i=1, imm1
                    wo(i, 1, k)=wo(i, 1, k)-cr(1, np1, k)*vb(i, np1, iv)
216             continue
                end do
            end do
        end do
        !
        !     case m = 1 through nlat-1
        !
        if (mmax < 2) goto 950
        do mp1=2, mmax
            m = mp1-1
            mp2 = mp1+1
            call sphere_aux%vbin(0, nlat, nlon, m, vb, iv, wvbin)
            call sphere_aux%wbin(0, nlat, nlon, m, wb, iw, wwbin)
            if (mp1 > ndo1) goto 226
            do k=1, nt
                do np1=mp1, ndo1, 2
                    do i=1, imm1
                        ve(i, 2*mp1-2, k) = ve(i, 2*mp1-2, k)-ci(mp1, np1, k)*wb(i, np1, iw)
                        ve(i, 2*mp1-1, k) = ve(i, 2*mp1-1, k)+cr(mp1, np1, k)*wb(i, np1, iw)
                        wo(i, 2*mp1-2, k) = wo(i, 2*mp1-2, k)-cr(mp1, np1, k)*vb(i, np1, iv)
                        wo(i, 2*mp1-1, k) = wo(i, 2*mp1-1, k)-ci(mp1, np1, k)*vb(i, np1, iv)
223                 continue
                    end do
                    if (mlat == 0) goto 224
                    ve(imid, 2*mp1-2, k) = ve(imid, 2*mp1-2, k) &
                        -ci(mp1, np1, k)*wb(imid, np1, iw)
                    ve(imid, 2*mp1-1, k) = ve(imid, 2*mp1-1, k) &
                        +cr(mp1, np1, k)*wb(imid, np1, iw)
224             continue
                end do
225         continue
            end do
226         if (mp2 > ndo2) goto 230
            do k=1, nt
                do np1=mp2, ndo2, 2
                    do i=1, imm1
                        vo(i, 2*mp1-2, k) = vo(i, 2*mp1-2, k)-ci(mp1, np1, k)*wb(i, np1, iw)
                        vo(i, 2*mp1-1, k) = vo(i, 2*mp1-1, k)+cr(mp1, np1, k)*wb(i, np1, iw)
                        we(i, 2*mp1-2, k) = we(i, 2*mp1-2, k)-cr(mp1, np1, k)*vb(i, np1, iv)
                        we(i, 2*mp1-1, k) = we(i, 2*mp1-1, k)-ci(mp1, np1, k)*vb(i, np1, iv)
227                 continue
                    end do
                    if (mlat == 0) goto 228
                    we(imid, 2*mp1-2, k) = we(imid, 2*mp1-2, k) &
                        -cr(mp1, np1, k)*vb(imid, np1, iv)
                    we(imid, 2*mp1-1, k) = we(imid, 2*mp1-1, k) &
                        -ci(mp1, np1, k)*vb(imid, np1, iv)
228             continue
                end do
229         continue
            end do
230     continue
        end do
        goto 950
        !
        !     case ityp=3   v even,  w odd
        !
300     call sphere_aux%vbin(0, nlat, nlon, 0, vb, iv, wvbin)
        !
        !     case m = 0
        !
        do k=1, nt
            do np1=2, ndo2, 2
                do i=1, imid
                    ve(i, 1, k)=ve(i, 1, k)+br(1, np1, k)*vb(i, np1, iv)
315             continue
                end do
            end do
        end do
        do k=1, nt
            do np1=3, ndo1, 2
                do i=1, imm1
                    wo(i, 1, k)=wo(i, 1, k)-cr(1, np1, k)*vb(i, np1, iv)
316             continue
                end do
            end do
        end do
        !
        !     case m = 1 through nlat-1
        !
        if (mmax < 2) goto 950
        do mp1=2, mmax
            m = mp1-1
            mp2 = mp1+1
            call sphere_aux%vbin(0, nlat, nlon, m, vb, iv, wvbin)
            call sphere_aux%wbin(0, nlat, nlon, m, wb, iw, wwbin)
            if (mp1 > ndo1) goto 326
            do k=1, nt
                do np1=mp1, ndo1, 2
                    do i=1, imm1
                        ve(i, 2*mp1-2, k) = ve(i, 2*mp1-2, k)-ci(mp1, np1, k)*wb(i, np1, iw)
                        ve(i, 2*mp1-1, k) = ve(i, 2*mp1-1, k)+cr(mp1, np1, k)*wb(i, np1, iw)
                        wo(i, 2*mp1-2, k) = wo(i, 2*mp1-2, k)-cr(mp1, np1, k)*vb(i, np1, iv)
                        wo(i, 2*mp1-1, k) = wo(i, 2*mp1-1, k)-ci(mp1, np1, k)*vb(i, np1, iv)
323                 continue
                    end do
                    if (mlat == 0) goto 324
                    ve(imid, 2*mp1-2, k) = ve(imid, 2*mp1-2, k) &
                        -ci(mp1, np1, k)*wb(imid, np1, iw)
                    ve(imid, 2*mp1-1, k) = ve(imid, 2*mp1-1, k) &
                        +cr(mp1, np1, k)*wb(imid, np1, iw)
324             continue
                end do
325         continue
            end do
326         if (mp2 > ndo2) goto 330
            do k=1, nt
                do np1=mp2, ndo2, 2
                    do i=1, imm1
                        ve(i, 2*mp1-2, k) = ve(i, 2*mp1-2, k)+br(mp1, np1, k)*vb(i, np1, iv)
                        ve(i, 2*mp1-1, k) = ve(i, 2*mp1-1, k)+bi(mp1, np1, k)*vb(i, np1, iv)
                        wo(i, 2*mp1-2, k) = wo(i, 2*mp1-2, k)-bi(mp1, np1, k)*wb(i, np1, iw)
                        wo(i, 2*mp1-1, k) = wo(i, 2*mp1-1, k)+br(mp1, np1, k)*wb(i, np1, iw)
327                 continue
                    end do
                    if (mlat == 0) goto 328
                    ve(imid, 2*mp1-2, k) = ve(imid, 2*mp1-2, k) &
                        +br(mp1, np1, k)*vb(imid, np1, iv)
                    ve(imid, 2*mp1-1, k) = ve(imid, 2*mp1-1, k) &
                        +bi(mp1, np1, k)*vb(imid, np1, iv)
328             continue
                end do
329         continue
            end do
330     continue
        end do
        goto 950
        !
        !     case ityp=4   v even,  w odd, and both cr and ci equal zero
        !
400     call sphere_aux%vbin(1, nlat, nlon, 0, vb, iv, wvbin)
        !
        !     case m = 0
        !
        do k=1, nt
            do np1=2, ndo2, 2
                do i=1, imid
                    ve(i, 1, k)=ve(i, 1, k)+br(1, np1, k)*vb(i, np1, iv)
415             continue
                end do
            end do
        end do
        !
        !     case m = 1 through nlat-1
        !
        if (mmax < 2) goto 950
        do mp1=2, mmax
            m = mp1-1
            mp2 = mp1+1
            call sphere_aux%vbin(1, nlat, nlon, m, vb, iv, wvbin)
            call sphere_aux%wbin(1, nlat, nlon, m, wb, iw, wwbin)
            if (mp2 > ndo2) goto 430
            do k=1, nt
                do np1=mp2, ndo2, 2
                    do i=1, imm1
                        ve(i, 2*mp1-2, k) = ve(i, 2*mp1-2, k)+br(mp1, np1, k)*vb(i, np1, iv)
                        ve(i, 2*mp1-1, k) = ve(i, 2*mp1-1, k)+bi(mp1, np1, k)*vb(i, np1, iv)
                        wo(i, 2*mp1-2, k) = wo(i, 2*mp1-2, k)-bi(mp1, np1, k)*wb(i, np1, iw)
                        wo(i, 2*mp1-1, k) = wo(i, 2*mp1-1, k)+br(mp1, np1, k)*wb(i, np1, iw)
427                 continue
                    end do
                    if (mlat == 0) goto 428
                    ve(imid, 2*mp1-2, k) = ve(imid, 2*mp1-2, k) &
                        +br(mp1, np1, k)*vb(imid, np1, iv)
                    ve(imid, 2*mp1-1, k) = ve(imid, 2*mp1-1, k) &
                        +bi(mp1, np1, k)*vb(imid, np1, iv)
428             continue
                end do
429         continue
            end do
430     continue
        end do
        goto 950
        !
        !     case ityp=5   v even,  w odd,     br and bi equal zero
        !
500     call sphere_aux%vbin(2, nlat, nlon, 0, vb, iv, wvbin)
        !
        !     case m = 0
        !
        do k=1, nt
            do np1=3, ndo1, 2
                do i=1, imm1
                    wo(i, 1, k)=wo(i, 1, k)-cr(1, np1, k)*vb(i, np1, iv)
516             continue
                end do
            end do
        end do
        !
        !     case m = 1 through nlat-1
        !
        if (mmax < 2) goto 950
        do mp1=2, mmax
            m = mp1-1
            mp2 = mp1+1
            call sphere_aux%vbin(2, nlat, nlon, m, vb, iv, wvbin)
            call sphere_aux%wbin(2, nlat, nlon, m, wb, iw, wwbin)
            if (mp1 > ndo1) goto 530
            do k=1, nt
                do np1=mp1, ndo1, 2
                    do i=1, imm1
                        ve(i, 2*mp1-2, k) = ve(i, 2*mp1-2, k)-ci(mp1, np1, k)*wb(i, np1, iw)
                        ve(i, 2*mp1-1, k) = ve(i, 2*mp1-1, k)+cr(mp1, np1, k)*wb(i, np1, iw)
                        wo(i, 2*mp1-2, k) = wo(i, 2*mp1-2, k)-cr(mp1, np1, k)*vb(i, np1, iv)
                        wo(i, 2*mp1-1, k) = wo(i, 2*mp1-1, k)-ci(mp1, np1, k)*vb(i, np1, iv)
523                 continue
                    end do
                    if (mlat == 0) goto 524
                    ve(imid, 2*mp1-2, k) = ve(imid, 2*mp1-2, k) &
                        -ci(mp1, np1, k)*wb(imid, np1, iw)
                    ve(imid, 2*mp1-1, k) = ve(imid, 2*mp1-1, k) &
                        +cr(mp1, np1, k)*wb(imid, np1, iw)
524             continue
                end do
525         continue
            end do
530     continue
        end do
        goto 950
        !
        !     case ityp=6   v odd  ,  w even
        !
600     call sphere_aux%vbin(0, nlat, nlon, 0, vb, iv, wvbin)
        !
        !     case m = 0
        !
        do k=1, nt
            do np1=2, ndo2, 2
                do i=1, imid
                    we(i, 1, k)=we(i, 1, k)-cr(1, np1, k)*vb(i, np1, iv)
615             continue
                end do
            end do
        end do
        do k=1, nt
            do np1=3, ndo1, 2
                do i=1, imm1
                    vo(i, 1, k)=vo(i, 1, k)+br(1, np1, k)*vb(i, np1, iv)
616             continue
                end do
            end do
        end do
        !
        !     case m = 1 through nlat-1
        !
        if (mmax < 2) goto 950
        do mp1=2, mmax
            m = mp1-1
            mp2 = mp1+1
            call sphere_aux%vbin(0, nlat, nlon, m, vb, iv, wvbin)
            call sphere_aux%wbin(0, nlat, nlon, m, wb, iw, wwbin)
            if (mp1 > ndo1) goto 626
            do k=1, nt
                do np1=mp1, ndo1, 2
                    do i=1, imm1
                        vo(i, 2*mp1-2, k) = vo(i, 2*mp1-2, k)+br(mp1, np1, k)*vb(i, np1, iv)
                        vo(i, 2*mp1-1, k) = vo(i, 2*mp1-1, k)+bi(mp1, np1, k)*vb(i, np1, iv)
                        we(i, 2*mp1-2, k) = we(i, 2*mp1-2, k)-bi(mp1, np1, k)*wb(i, np1, iw)
                        we(i, 2*mp1-1, k) = we(i, 2*mp1-1, k)+br(mp1, np1, k)*wb(i, np1, iw)
623                 continue
                    end do
                    if (mlat == 0) goto 624
                    we(imid, 2*mp1-2, k) = we(imid, 2*mp1-2, k) &
                        -bi(mp1, np1, k)*wb(imid, np1, iw)
                    we(imid, 2*mp1-1, k) = we(imid, 2*mp1-1, k) &
                        +br(mp1, np1, k)*wb(imid, np1, iw)
624             continue
                end do
625         continue
            end do
626         if (mp2 > ndo2) goto 630
            do k=1, nt
                do np1=mp2, ndo2, 2
                    do i=1, imm1
                        vo(i, 2*mp1-2, k) = vo(i, 2*mp1-2, k)-ci(mp1, np1, k)*wb(i, np1, iw)
                        vo(i, 2*mp1-1, k) = vo(i, 2*mp1-1, k)+cr(mp1, np1, k)*wb(i, np1, iw)
                        we(i, 2*mp1-2, k) = we(i, 2*mp1-2, k)-cr(mp1, np1, k)*vb(i, np1, iv)
                        we(i, 2*mp1-1, k) = we(i, 2*mp1-1, k)-ci(mp1, np1, k)*vb(i, np1, iv)
627                 continue
                    end do
                    if (mlat == 0) goto 628
                    we(imid, 2*mp1-2, k) = we(imid, 2*mp1-2, k) &
                        -cr(mp1, np1, k)*vb(imid, np1, iv)
                    we(imid, 2*mp1-1, k) = we(imid, 2*mp1-1, k) &
                        -ci(mp1, np1, k)*vb(imid, np1, iv)
628             continue
                end do
629         continue
            end do
630     continue
        end do
        goto 950
        !
        !     case ityp=7   v odd, w even   cr and ci equal zero
        !
700     call sphere_aux%vbin(2, nlat, nlon, 0, vb, iv, wvbin)
        !
        !     case m = 0
        !
        do k=1, nt
            do np1=3, ndo1, 2
                do i=1, imm1
                    vo(i, 1, k)=vo(i, 1, k)+br(1, np1, k)*vb(i, np1, iv)
716             continue
                end do
            end do
        end do
        !
        !     case m = 1 through nlat-1
        !
        if (mmax < 2) goto 950
        do mp1=2, mmax
            m = mp1-1
            mp2 = mp1+1
            call sphere_aux%vbin(2, nlat, nlon, m, vb, iv, wvbin)
            call sphere_aux%wbin(2, nlat, nlon, m, wb, iw, wwbin)
            if (mp1 > ndo1) goto 730
            do k=1, nt
                do np1=mp1, ndo1, 2
                    do i=1, imm1
                        vo(i, 2*mp1-2, k) = vo(i, 2*mp1-2, k)+br(mp1, np1, k)*vb(i, np1, iv)
                        vo(i, 2*mp1-1, k) = vo(i, 2*mp1-1, k)+bi(mp1, np1, k)*vb(i, np1, iv)
                        we(i, 2*mp1-2, k) = we(i, 2*mp1-2, k)-bi(mp1, np1, k)*wb(i, np1, iw)
                        we(i, 2*mp1-1, k) = we(i, 2*mp1-1, k)+br(mp1, np1, k)*wb(i, np1, iw)
723                 continue
                    end do
                    if (mlat == 0) goto 724
                    we(imid, 2*mp1-2, k) = we(imid, 2*mp1-2, k) &
                        -bi(mp1, np1, k)*wb(imid, np1, iw)
                    we(imid, 2*mp1-1, k) = we(imid, 2*mp1-1, k) &
                        +br(mp1, np1, k)*wb(imid, np1, iw)
724             continue
                end do
725         continue
            end do
730     continue
        end do
        goto 950
        !
        !     case ityp=8   v odd,  w even   br and bi equal zero
        !
800     call sphere_aux%vbin(1, nlat, nlon, 0, vb, iv, wvbin)
        !
        !     case m = 0
        !
        do k=1, nt
            do np1=2, ndo2, 2
                do i=1, imid
                    we(i, 1, k)=we(i, 1, k)-cr(1, np1, k)*vb(i, np1, iv)
815             continue
                end do
            end do
        end do
        !
        !     case m = 1 through nlat-1
        !
        if (mmax < 2) goto 950
        do mp1=2, mmax
            m = mp1-1
            mp2 = mp1+1
            call sphere_aux%vbin(1, nlat, nlon, m, vb, iv, wvbin)
            call sphere_aux%wbin(1, nlat, nlon, m, wb, iw, wwbin)
            if (mp2 > ndo2) goto 830
            do k=1, nt
                do np1=mp2, ndo2, 2
                    do i=1, imm1
                        vo(i, 2*mp1-2, k) = vo(i, 2*mp1-2, k)-ci(mp1, np1, k)*wb(i, np1, iw)
                        vo(i, 2*mp1-1, k) = vo(i, 2*mp1-1, k)+cr(mp1, np1, k)*wb(i, np1, iw)
                        we(i, 2*mp1-2, k) = we(i, 2*mp1-2, k)-cr(mp1, np1, k)*vb(i, np1, iv)
                        we(i, 2*mp1-1, k) = we(i, 2*mp1-1, k)-ci(mp1, np1, k)*vb(i, np1, iv)
827                 continue
                    end do
                    if (mlat == 0) goto 828
                    we(imid, 2*mp1-2, k) = we(imid, 2*mp1-2, k) &
                        -cr(mp1, np1, k)*vb(imid, np1, iv)
                    we(imid, 2*mp1-1, k) = we(imid, 2*mp1-1, k) &
                        -ci(mp1, np1, k)*vb(imid, np1, iv)
828             continue
                end do
829         continue
            end do
830     continue
        end do
        950 do k=1, nt
            call sphere_aux%hfft%backward(idv, nlon, ve(1, 1, k), idv, wrfft, vb)
            call sphere_aux%hfft%backward(idv, nlon, we(1, 1, k), idv, wrfft, vb)
14      continue
        end do
        if (ityp > 2) goto 12
        do k=1, nt
            do j=1, nlon
                do i=1, imm1
                    v(i, j, k) = HALF *(ve(i, j, k)+vo(i, j, k))
                    w(i, j, k) = HALF *(we(i, j, k)+wo(i, j, k))
                    v(nlp1-i, j, k) = HALF *(ve(i, j, k)-vo(i, j, k))
                    w(nlp1-i, j, k) = HALF *(we(i, j, k)-wo(i, j, k))
60              continue
                end do
            end do
        end do
        goto 13
        12 do k=1, nt
            do j=1, nlon
                do i=1, imm1
                    v(i, j, k) = HALF * ve(i, j, k)
                    w(i, j, k) = HALF * we(i, j, k)
11              continue
                end do
            end do
        end do
13      if (mlat == 0) return
        do k=1, nt
            do j=1, nlon
                v(imid, j, k) = HALF * ve(imid, j, k)
                w(imid, j, k) = HALF * we(imid, j, k)
65          continue
            end do
        end do

    end subroutine vhsec_lower_routine

end submodule vector_synthesis_regular_grid
