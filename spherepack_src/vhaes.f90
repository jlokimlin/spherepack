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
! ... file vhaes.f
!
!     this file contains code and documentation for subroutines
!     vhaes and vhaesi
!
! ... files which must be loaded with vhaes.f
!
!     sphcom.f, hrfft.f
!
!                                                                              
!     subroutine vhaes(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci,
!    +                 mdab, ndab, wvhaes, lvhaes, work, lwork, ierror)
!
!     subroutine vhaes performs the vector spherical harmonic analysis
!     on the vector field (v, w) and stores the result in the arrays
!     br, bi, cr, and ci. v(i, j) and w(i, j) are the colatitudinal
!     (measured from the north pole) and east longitudinal components
!     respectively, located at colatitude theta(i) = (i-1)*pi/(nlat-1)
!     and longitude phi(j) = (j-1)*2*pi/nlon. the spectral
!     representation of (v, w) is given at output parameters v, w in
!     subroutine vhses.  
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
!     ityp   = 0  no symmetries exist about the equator. the analysis
!                 is performed on the entire sphere.  i.e. on the
!                 arrays v(i, j), w(i, j) for i=1, ..., nlat and
!                 j=1, ..., nlon.
!
!            = 1  no symmetries exist about the equator. the analysis
!                 is performed on the entire sphere.  i.e. on the
!                 arrays v(i, j), w(i, j) for i=1, ..., nlat and
!                 j=1, ..., nlon. the curl of (v, w) is zero. that is,
!                 (d/dtheta (sin(theta) w) - dv/dphi)/sin(theta) = 0. 
!                 the coefficients cr and ci are zero.
!
!            = 2  no symmetries exist about the equator. the analysis
!                 is performed on the entire sphere.  i.e. on the
!                 arrays v(i, j), w(i, j) for i=1, ..., nlat and
!                 j=1, ..., nlon. the divergence of (v, w) is zero. i.e.,
!                 (d/dtheta (sin(theta) v) + dw/dphi)/sin(theta) = 0. 
!                 the coefficients br and bi are zero.
!
!            = 3  v is symmetric and w is antisymmetric about the 
!                 equator. the analysis is performed on the northern
!                 hemisphere only.  i.e., if nlat is odd the analysis
!                 is performed on the arrays v(i, j), w(i, j) for
!                 i=1, ..., (nlat+1)/2 and j=1, ..., nlon. if nlat is
!                 even the analysis is performed on the the arrays
!                 v(i, j), w(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
!
!            = 4  v is symmetric and w is antisymmetric about the 
!                 equator. the analysis is performed on the northern
!                 hemisphere only.  i.e., if nlat is odd the analysis
!                 is performed on the arrays v(i, j), w(i, j) for
!                 i=1, ..., (nlat+1)/2 and j=1, ..., nlon. if nlat is
!                 even the analysis is performed on the the arrays
!                 v(i, j), w(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
!                 the curl of (v, w) is zero. that is,
!                 (d/dtheta (sin(theta) w) - dv/dphi)/sin(theta) = 0. 
!                 the coefficients cr and ci are zero.
!
!            = 5  v is symmetric and w is antisymmetric about the 
!                 equator. the analysis is performed on the northern
!                 hemisphere only.  i.e., if nlat is odd the analysis
!                 is performed on the arrays v(i, j), w(i, j) for
!                 i=1, ..., (nlat+1)/2 and j=1, ..., nlon. if nlat is
!                 even the analysis is performed on the the arrays
!                 v(i, j), w(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
!                 the divergence of (v, w) is zero. i.e.,
!                 (d/dtheta (sin(theta) v) + dw/dphi)/sin(theta) = 0. 
!                 the coefficients br and bi are zero.
!
!            = 6  v is antisymmetric and w is symmetric about the 
!                 equator. the analysis is performed on the northern
!                 hemisphere only.  i.e., if nlat is odd the analysis
!                 is performed on the arrays v(i, j), w(i, j) for
!                 i=1, ..., (nlat+1)/2 and j=1, ..., nlon. if nlat is
!                 even the analysis is performed on the the arrays
!                 v(i, j), w(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
!
!            = 7  v is antisymmetric and w is symmetric about the 
!                 equator. the analysis is performed on the northern
!                 hemisphere only.  i.e., if nlat is odd the analysis
!                 is performed on the arrays v(i, j), w(i, j) for
!                 i=1, ..., (nlat+1)/2 and j=1, ..., nlon. if nlat is
!                 even the analysis is performed on the the arrays
!                 v(i, j), w(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
!                 the curl of (v, w) is zero. that is,
!                 (d/dtheta (sin(theta) w) - dv/dphi)/sin(theta) = 0. 
!                 the coefficients cr and ci are zero.
!
!            = 8  v is antisymmetric and w is symmetric about the 
!                 equator. the analysis is performed on the northern
!                 hemisphere only.  i.e., if nlat is odd the analysis
!                 is performed on the arrays v(i, j), w(i, j) for
!                 i=1, ..., (nlat+1)/2 and j=1, ..., nlon. if nlat is
!                 even the analysis is performed on the the arrays
!                 v(i, j), w(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
!                 the divergence of (v, w) is zero. i.e.,
!                 (d/dtheta (sin(theta) v) + dw/dphi)/sin(theta) = 0. 
!                 the coefficients br and bi are zero.
!
!
!     nt     the number of analyses.  in the program that calls vhaes,
!            the arrays v, w, br, bi, cr, and ci can be three dimensional
!            in which case multiple analyses will be performed.
!            the third index is the analysis index which assumes the 
!            values k=1, ..., nt.  for a single analysis set nt=1. the
!            discription of the remaining parameters is simplified
!            by assuming that nt=1 or that all the arrays are two
!            dimensional.
!
!     v, w    two or three dimensional arrays (see input parameter nt)
!            that contain the vector function to be analyzed.
!            v is the colatitudnal component and w is the east 
!            longitudinal component. v(i, j), w(i, j) contain the
!            components at colatitude theta(i) = (i-1)*pi/(nlat-1)
!            and longitude phi(j) = (j-1)*2*pi/nlon. the index ranges
!            are defined above at the input parameter ityp.
!
!     idvw   the first dimension of the arrays v, w as it appears in
!            the program that calls vhaes. if ityp .le. 2 then idvw
!            must be at least nlat.  if ityp .gt. 2 and nlat is
!            even then idvw must be at least nlat/2. if ityp .gt. 2
!            and nlat is odd then idvw must be at least (nlat+1)/2.
!
!     jdvw   the second dimension of the arrays v, w as it appears in
!            the program that calls vhaes. jdvw must be at least nlon.
!
!     mdab   the first dimension of the arrays br, bi, cr, and ci as it
!            appears in the program that calls vhaes. mdab must be at
!            least min(nlat, nlon/2) if nlon is even or at least
!            min(nlat, (nlon+1)/2) if nlon is odd.
!
!     ndab   the second dimension of the arrays br, bi, cr, and ci as it
!            appears in the program that calls vhaes. ndab must be at
!            least nlat.
!
!     lvhaes an array which must be initialized by subroutine vhaesi.
!            once initialized, wvhaes can be used repeatedly by vhaes
!            as long as nlon and nlat remain unchanged.  wvhaes must
!            not be altered between calls of vhaes.
!
!     lvhaes the dimension of the array wvhaes as it appears in the
!            program that calls vhaes. define
!
!               l1 = min(nlat, nlon/2) if nlon is even or
!               l1 = min(nlat, (nlon+1)/2) if nlon is odd
!
!            and
!
!               l2 = nlat/2        if nlat is even or
!               l2 = (nlat+1)/2    if nlat is odd
!
!            then lvhaes must be at least
!
!            l1*l2(nlat+nlat-l1+1)+nlon+15
!
!
!     work   a work array that does not have to be saved.
!
!     lwork  the dimension of the array work as it appears in the
!            program that calls vhaes. define
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
!     br, bi  two or three dimensional arrays (see input parameter nt)
!     cr, ci  that contain the vector spherical harmonic coefficients
!            in the spectral representation of v(i, j) and w(i, j) given
!            in the discription of subroutine vhses. br(mp1, np1),
!            bi(mp1, np1), cr(mp1, np1), and ci(mp1, np1) are computed
!            for mp1=1, ..., mmax and np1=mp1, ..., nlat except for np1=nlat
!            and odd mp1. mmax=min(nlat, nlon/2) if nlon is even or
!            mmax=min(nlat, (nlon+1)/2) if nlon is odd.
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
!            = 9  error in the specification of lvhaes
!            = 10 error in the specification of lwork
!
! ********************************************************
!
!     subroutine vhaesi(nlat, nlon, wvhaes, lvhaes, work, lwork, dwork,
!    +                  ldwork, ierror)
!
!     subroutine vhaesi initializes the array wvhaes which can then be
!     used repeatedly by subroutine vhaes until nlat or nlon is changed.
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
!     lvhaes the dimension of the array wvhaes as it appears in the
!            program that calls vhaes. define
!
!               l1 = min(nlat, nlon/2) if nlon is even or
!               l1 = min(nlat, (nlon+1)/2) if nlon is odd
!
!            and
!
!               l2 = nlat/2        if nlat is even or
!               l2 = (nlat+1)/2    if nlat is odd
!
!            then lvhaes must be at least
!
!               l1*l2*(nlat+nlat-l1+1)+nlon+15
!
!
!     work   a work array that does not have to be saved.
!
!     lwork  the dimension of the array work as it appears in the
!            program that calls vhaes. lwork must be at least
!
!              3*(max(l1-2, 0)*(nlat+nlat-l1-1))/2+5*l2*nlat
!
!     dwork  an unsaved real work space
!
!     ldwork the length of the array dwork as it appears in the
!            program that calls vhaesi.  ldwork must be at least
!            2*(nlat+1)
!
!
!     **************************************************************
!
!     output parameters
!
!     wvhaes an array which is initialized for use by subroutine vhaes.
!            once initialized, wvhaes can be used repeatedly by vhaes
!            as long as nlat or nlon remain unchanged.  wvhaes must not
!            be altered between calls of vhaes.
!
!
!     ierror = 0  no errors
!            = 1  error in the specification of nlat
!            = 2  error in the specification of nlon
!            = 3  error in the specification of lvhaes
!            = 4  error in the specification of lwork
!            = 5  error in the specification of ldwork
!
!
subroutine vhaes(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, &
    mdab, ndab, wvhaes, lvhaes, work, lwork, ierror)
    dimension v(idvw, jdvw, 1), w(idvw, jdvw, 1), br(mdab, ndab, 1), &
        bi(mdab, ndab, 1), cr(mdab, ndab, 1), ci(mdab, ndab, 1), &
        work(1), wvhaes(1)
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
    if (lvhaes < lzimn+lzimn+nlon+15) return
    ierror = 10
    idv = nlat
    if (ityp > 2) idv = imid
    lnl = nt*idv*nlon
    if (lwork < lnl+lnl+idv*nlon) return
    ierror = 0
    ist = 0
    if (ityp <= 2) ist = imid
    iw1 = ist+1
    iw2 = lnl+1
    iw3 = iw2+ist
    iw4 = iw2+lnl
    jw1 = lzimn+1
    jw2 = jw1+lzimn

    call vhaes1(nlat, nlon, ityp, nt, imid, idvw, jdvw, v, w, mdab, ndab, &
        br, bi, cr, ci, idv, work, work(iw1), work(iw2), work(iw3), &
        work(iw4), idz, wvhaes, wvhaes(jw1), wvhaes(jw2))

end subroutine vhaes



subroutine vhaes1(nlat, nlon, ityp, nt, imid, idvw, jdvw, v, w, mdab, &
    ndab, br, bi, cr, ci, idv, ve, vo, we, wo, work, idz, zv, zw, wrfft)
    dimension v(idvw, jdvw, 1), w(idvw, jdvw, 1), br(mdab, ndab, 1), &
        bi(mdab, ndab, 1), cr(mdab, ndab, 1), ci(mdab, ndab, 1), &
        ve(idv, nlon, 1), vo(idv, nlon, 1), we(idv, nlon, 1), &
        wo(idv, nlon, 1), work(1), wrfft(1), &
        zv(idz, 1), zw(idz, 1)

    nlp1 = nlat+1
    tsn = 2.0/nlon
    fsn = 4.0/nlon
    mlat = mod(nlat, 2)
    mlon = mod(nlon, 2)
    mmax = min(nlat, (nlon+1)/2)
    imm1 = imid

    select case (mlat)
        case (0)
            imm1 = imid
            ndo1 = nlat
            ndo2 = nlat-1
        case default
            imm1 = imid-1
            ndo1 = nlat-1
            ndo2 = nlat
    end select

    if (ityp > 2) then
        go to 3
    end if

    do k=1, nt
        do i=1, imm1
            do j=1, nlon
                ve(i, j, k) = tsn*(v(i, j, k)+v(nlp1-i, j, k))
                vo(i, j, k) = tsn*(v(i, j, k)-v(nlp1-i, j, k))
                we(i, j, k) = tsn*(w(i, j, k)+w(nlp1-i, j, k))
                wo(i, j, k) = tsn*(w(i, j, k)-w(nlp1-i, j, k))
            end do
        end do
    end do

    go to 2

    3 do k=1, nt
        do i=1, imm1
            do j=1, nlon
                ve(i, j, k) = fsn*v(i, j, k)
                vo(i, j, k) = fsn*v(i, j, k)
                we(i, j, k) = fsn*w(i, j, k)
                wo(i, j, k) = fsn*w(i, j, k)
            end do
        end do
    end do

2   if (mlat == 0) then
        go to 7
    end if

    do k=1, nt
        do j=1, nlon
            ve(imid, j, k) = tsn*v(imid, j, k)
            we(imid, j, k) = tsn*w(imid, j, k)
        end do
    end do

    7 do k=1, nt
        call hrfftf(idv, nlon, ve(1, 1, k), idv, wrfft, work)
        call hrfftf(idv, nlon, we(1, 1, k), idv, wrfft, work)
    end do

    !
    !==> Set polar coefficients to zero
    !
    select case (ityp)
        case (0,1,3,4,6,7)
            do k=1, nt
                do mp1=1, mmax
                    do np1=mp1, nlat
                        br(mp1, np1, k) = 0.0
                        bi(mp1, np1, k) = 0.0
                    end do
                end do
            end do
    end select

    !
    !==> Set azimuthal coefficients to zero
    !
    select case (ityp)
        case (0,2,3,5,6,8)
            do k=1, nt
                do mp1=1, mmax
                    do np1=mp1, nlat
                        cr(mp1, np1, k) = 0.0
                        ci(mp1, np1, k) = 0.0
                    end do
                end do
            end do
    end select

    select case (ityp)
        case (0)
            !
            !==> case ityp=0,  no symmetries
            !
            ! case m=0
            !
            1 do k=1, nt
                do i=1, imid
                    do np1=2, ndo2, 2
                        br(1, np1, k) = br(1, np1, k)+zv(np1, i)*ve(i, 1, k)
                        cr(1, np1, k) = cr(1, np1, k)-zv(np1, i)*we(i, 1, k)
                    end do
                end do
            end do

            do k=1, nt
                do i=1, imm1
                    do np1=3, ndo1, 2
                        br(1, np1, k) = br(1, np1, k)+zv(np1, i)*vo(i, 1, k)
                        cr(1, np1, k) = cr(1, np1, k)-zv(np1, i)*wo(i, 1, k)
                    end do
                end do
            end do
            !
            !     case m = 1 through nlat-1
            !
            if (mmax < 2) then
                return
            end if

            do mp1=2, mmax
                m = mp1-1
                mb = m*(nlat-1)-(m*(m-1))/2
                mp2 = mp1+1

                if (mp1 > ndo1) then
                    go to 17
                end if

                do k=1, nt
                    do i=1, imm1
                        do np1=mp1, ndo1, 2
                            br(mp1, np1, k) = br(mp1, np1, k)+zv(np1+mb, i)*vo(i, 2*mp1-2, k) &
                                +zw(np1+mb, i)*we(i, 2*mp1-1, k)
                            bi(mp1, np1, k) = bi(mp1, np1, k)+zv(np1+mb, i)*vo(i, 2*mp1-1, k) &
                                -zw(np1+mb, i)*we(i, 2*mp1-2, k)
                            cr(mp1, np1, k) = cr(mp1, np1, k)-zv(np1+mb, i)*wo(i, 2*mp1-2, k) &
                                +zw(np1+mb, i)*ve(i, 2*mp1-1, k)
                            ci(mp1, np1, k) = ci(mp1, np1, k)-zv(np1+mb, i)*wo(i, 2*mp1-1, k) &
                                -zw(np1+mb, i)*ve(i, 2*mp1-2, k)
                        end do
                    end do
                end do

                if (mlat == 0) then
                    go to 17
                end if

                do k=1, nt
                    do np1=mp1, ndo1, 2
                        br(mp1, np1, k) = br(mp1, np1, k)+zw(np1+mb, imid)*we(imid, 2*mp1-1, k)
                        bi(mp1, np1, k) = bi(mp1, np1, k)-zw(np1+mb, imid)*we(imid, 2*mp1-2, k)
                        cr(mp1, np1, k) = cr(mp1, np1, k)+zw(np1+mb, imid)*ve(imid, 2*mp1-1, k)
                        ci(mp1, np1, k) = ci(mp1, np1, k)-zw(np1+mb, imid)*ve(imid, 2*mp1-2, k)
                    end do
                end do

17              if (mp2 > ndo2) then
                    exit
                end if

                do k=1, nt
                    do i=1, imm1
                        do np1=mp2, ndo2, 2
                            br(mp1, np1, k) = br(mp1, np1, k)+zv(np1+mb, i)*ve(i, 2*mp1-2, k) &
                                +zw(np1+mb, i)*wo(i, 2*mp1-1, k)
                            bi(mp1, np1, k) = bi(mp1, np1, k)+zv(np1+mb, i)*ve(i, 2*mp1-1, k) &
                                -zw(np1+mb, i)*wo(i, 2*mp1-2, k)
                            cr(mp1, np1, k) = cr(mp1, np1, k)-zv(np1+mb, i)*we(i, 2*mp1-2, k) &
                                +zw(np1+mb, i)*vo(i, 2*mp1-1, k)
                            ci(mp1, np1, k) = ci(mp1, np1, k)-zv(np1+mb, i)*we(i, 2*mp1-1, k) &
                                -zw(np1+mb, i)*vo(i, 2*mp1-2, k)
                        end do
                    end do
                end do

                if (mlat == 0) then
                    exit
                end if

                do k=1, nt
                    do np1=mp2, ndo2, 2
                        br(mp1, np1, k) = br(mp1, np1, k)+zv(np1+mb, imid)*ve(imid, 2*mp1-2, k)
                        bi(mp1, np1, k) = bi(mp1, np1, k)+zv(np1+mb, imid)*ve(imid, 2*mp1-1, k)
                        cr(mp1, np1, k) = cr(mp1, np1, k)-zv(np1+mb, imid)*we(imid, 2*mp1-2, k)
                        ci(mp1, np1, k) = ci(mp1, np1, k)-zv(np1+mb, imid)*we(imid, 2*mp1-1, k)
                    end do
                end do
            end do
            return
        case (1)
            !
            !==> case ityp=1 ,  no symmetries but cr and ci equal zero
            !
            !    case m=0
            !
            do k=1, nt
                do i=1, imid
                    do np1=2, ndo2, 2
                        br(1, np1, k) = br(1, np1, k)+zv(np1, i)*ve(i, 1, k)
                    end do
                end do
            end do

            do k=1, nt
                do i=1, imm1
                    do np1=3, ndo1, 2
                        br(1, np1, k) = br(1, np1, k)+zv(np1, i)*vo(i, 1, k)
                    end do
                end do
            end do
            !
            !==> case m = 1 through nlat-1
            !
            if (mmax < 2) then
                return
            end if

            do mp1=2, mmax
                m = mp1-1
                mb = m*(nlat-1)-(m*(m-1))/2
                mp2 = mp1+1

                if (mp1 > ndo1) then
                    go to 117
                end if

                do k=1, nt
                    do i=1, imm1
                        do np1=mp1, ndo1, 2
                            br(mp1, np1, k) = br(mp1, np1, k)+zv(np1+mb, i)*vo(i, 2*mp1-2, k) &
                                +zw(np1+mb, i)*we(i, 2*mp1-1, k)
                            bi(mp1, np1, k) = bi(mp1, np1, k)+zv(np1+mb, i)*vo(i, 2*mp1-1, k) &
                                -zw(np1+mb, i)*we(i, 2*mp1-2, k)
                        end do
                    end do
                end do

                if (mlat == 0) then
                    go to 117
                end if

                do k=1, nt
                    do np1=mp1, ndo1, 2
                        br(mp1, np1, k) = br(mp1, np1, k)+zw(np1+mb, imid)*we(imid, 2*mp1-1, k)
                        bi(mp1, np1, k) = bi(mp1, np1, k)-zw(np1+mb, imid)*we(imid, 2*mp1-2, k)
                    end do
                end do

117             if (mp2 > ndo2) then
                    exit
                end if

                do k=1, nt
                    do i=1, imm1
                        do np1=mp2, ndo2, 2
                            br(mp1, np1, k) = br(mp1, np1, k)+zv(np1+mb, i)*ve(i, 2*mp1-2, k) &
                                +zw(np1+mb, i)*wo(i, 2*mp1-1, k)
                            bi(mp1, np1, k) = bi(mp1, np1, k)+zv(np1+mb, i)*ve(i, 2*mp1-1, k) &
                                -zw(np1+mb, i)*wo(i, 2*mp1-2, k)
                        end do
                    end do
                end do

                if (mlat == 0) then
                    exit
                end if

                do k=1, nt
                    do np1=mp2, ndo2, 2
                        br(mp1, np1, k) = br(mp1, np1, k)+zv(np1+mb, imid)*ve(imid, 2*mp1-2, k)
                        bi(mp1, np1, k) = bi(mp1, np1, k)+zv(np1+mb, imid)*ve(imid, 2*mp1-1, k)
                    end do
                end do
            end do
            return
        case (2)
            !
            !==> case ityp=2 ,  no symmetries but br and bi equal zero
            !
            !    case m=0
            !
            do k=1, nt
                do i=1, imid
                    do np1=2, ndo2, 2
                        cr(1, np1, k) = cr(1, np1, k)-zv(np1, i)*we(i, 1, k)
                    end do
                end do
            end do

            do k=1, nt
                do i=1, imm1
                    do np1=3, ndo1, 2
                        cr(1, np1, k) = cr(1, np1, k)-zv(np1, i)*wo(i, 1, k)
                    end do
                end do
            end do
            !
            !==> case m = 1 through nlat-1
            !
            if (mmax < 2) then
                return
            end if

            do mp1=2, mmax
                m = mp1-1
                mb = m*(nlat-1)-(m*(m-1))/2
                mp2 = mp1+1

                if (mp1 > ndo1) then
                    go to 217
                end if

                do k=1, nt
                    do i=1, imm1
                        do np1=mp1, ndo1, 2
                            cr(mp1, np1, k) = cr(mp1, np1, k)-zv(np1+mb, i)*wo(i, 2*mp1-2, k) &
                                +zw(np1+mb, i)*ve(i, 2*mp1-1, k)
                            ci(mp1, np1, k) = ci(mp1, np1, k)-zv(np1+mb, i)*wo(i, 2*mp1-1, k) &
                                -zw(np1+mb, i)*ve(i, 2*mp1-2, k)
                        end do
                    end do
                end do

                if (mlat == 0) then
                    go to 217
                end if

                do k=1, nt
                    do np1=mp1, ndo1, 2
                        cr(mp1, np1, k) = cr(mp1, np1, k)+zw(np1+mb, imid)*ve(imid, 2*mp1-1, k)
                        ci(mp1, np1, k) = ci(mp1, np1, k)-zw(np1+mb, imid)*ve(imid, 2*mp1-2, k)
                    end do
                end do

217             if (mp2 > ndo2) then
                    exit
                end if

                do k=1, nt
                    do i=1, imm1
                        do np1=mp2, ndo2, 2
                            cr(mp1, np1, k) = cr(mp1, np1, k)-zv(np1+mb, i)*we(i, 2*mp1-2, k) &
                                +zw(np1+mb, i)*vo(i, 2*mp1-1, k)
                            ci(mp1, np1, k) = ci(mp1, np1, k)-zv(np1+mb, i)*we(i, 2*mp1-1, k) &
                                -zw(np1+mb, i)*vo(i, 2*mp1-2, k)
                        end do
                    end do
                end do

                if (mlat == 0) then
                    exit
                end if

                do k=1, nt
                    do np1=mp2, ndo2, 2
                        cr(mp1, np1, k) = cr(mp1, np1, k)-zv(np1+mb, imid)*we(imid, 2*mp1-2, k)
                        ci(mp1, np1, k) = ci(mp1, np1, k)-zv(np1+mb, imid)*we(imid, 2*mp1-1, k)
                    end do
                end do
            end do
            return
        case (3)
            !
            !==> case ityp=3 ,  v even , w odd
            !
            !    case m=0
            !
            do k=1, nt
                do i=1, imid
                    do np1=2, ndo2, 2
                        br(1, np1, k) = br(1, np1, k)+zv(np1, i)*ve(i, 1, k)
                    end do
                end do
            end do

            do k=1, nt
                do i=1, imm1
                    do np1=3, ndo1, 2
                        cr(1, np1, k) = cr(1, np1, k)-zv(np1, i)*wo(i, 1, k)
                    end do
                end do
            end do
            !
            !==> case m = 1 through nlat-1
            !
            if (mmax < 2) then
                return
            end if

            do mp1=2, mmax
                m = mp1-1
                mb = m*(nlat-1)-(m*(m-1))/2
                mp2 = mp1+1

                if (mp1 > ndo1) then
                    go to 317
                end if

                do k=1, nt
                    do i=1, imm1
                        do np1=mp1, ndo1, 2
                            cr(mp1, np1, k) = cr(mp1, np1, k)-zv(np1+mb, i)*wo(i, 2*mp1-2, k) &
                                +zw(np1+mb, i)*ve(i, 2*mp1-1, k)
                            ci(mp1, np1, k) = ci(mp1, np1, k)-zv(np1+mb, i)*wo(i, 2*mp1-1, k) &
                                -zw(np1+mb, i)*ve(i, 2*mp1-2, k)
                        end do
                    end do
                end do

                if (mlat == 0) then
                    go to 317
                end if

                do k=1, nt
                    do np1=mp1, ndo1, 2
                        cr(mp1, np1, k) = cr(mp1, np1, k)+zw(np1+mb, imid)*ve(imid, 2*mp1-1, k)
                        ci(mp1, np1, k) = ci(mp1, np1, k)-zw(np1+mb, imid)*ve(imid, 2*mp1-2, k)
                    end do
                end do

317             if (mp2 > ndo2) then
                    exit
                end if

                do k=1, nt
                    do i=1, imm1
                        do np1=mp2, ndo2, 2
                            br(mp1, np1, k) = br(mp1, np1, k)+zv(np1+mb, i)*ve(i, 2*mp1-2, k) &
                                +zw(np1+mb, i)*wo(i, 2*mp1-1, k)
                            bi(mp1, np1, k) = bi(mp1, np1, k)+zv(np1+mb, i)*ve(i, 2*mp1-1, k) &
                                -zw(np1+mb, i)*wo(i, 2*mp1-2, k)
                        end do
                    end do
                end do

                if (mlat == 0) then
                    exit
                end if

                do  k=1, nt
                    do  np1=mp2, ndo2, 2
                        br(mp1, np1, k) = br(mp1, np1, k)+zv(np1+mb, imid)*ve(imid, 2*mp1-2, k)
                        bi(mp1, np1, k) = bi(mp1, np1, k)+zv(np1+mb, imid)*ve(imid, 2*mp1-1, k)
                    end do
                end do
            end do
            return
        case (4)
            !
            !==> case ityp=4 ,  v even, w odd, and cr and ci equal 0.
            !
            !    case m=0
            !
            do k=1, nt
                do i=1, imid
                    do np1=2, ndo2, 2
                        br(1, np1, k) = br(1, np1, k)+zv(np1, i)*ve(i, 1, k)
                    end do
                end do
            end do
            !
            !==> case m = 1 through nlat-1
            !
            if (mmax < 2) then
                return
            end if

            do mp1=2, mmax
                m = mp1-1
                mb = m*(nlat-1)-(m*(m-1))/2
                mp2 = mp1+1

                if (mp2 > ndo2) then
                    exit
                end if

                do k=1, nt
                    do i=1, imm1
                        do np1=mp2, ndo2, 2
                            br(mp1, np1, k) = br(mp1, np1, k)+zv(np1+mb, i)*ve(i, 2*mp1-2, k) &
                                +zw(np1+mb, i)*wo(i, 2*mp1-1, k)
                            bi(mp1, np1, k) = bi(mp1, np1, k)+zv(np1+mb, i)*ve(i, 2*mp1-1, k) &
                                -zw(np1+mb, i)*wo(i, 2*mp1-2, k)
                        end do
                    end do
                end do

                if (mlat == 0) then
                    exit
                end if

                do k=1, nt
                    do np1=mp2, ndo2, 2
                        br(mp1, np1, k) = br(mp1, np1, k)+zv(np1+mb, imid)*ve(imid, 2*mp1-2, k)
                        bi(mp1, np1, k) = bi(mp1, np1, k)+zv(np1+mb, imid)*ve(imid, 2*mp1-1, k)
                    end do
                end do
            end do
            return
        case (5)
            !
            !==> case ityp=5   v even, w odd, and br and bi equal zero
            !
            !    case m=0
            !
            do k=1, nt
                do i=1, imm1
                    do np1=3, ndo1, 2
                        cr(1, np1, k) = cr(1, np1, k)-zv(np1, i)*wo(i, 1, k)
                    end do
                end do
            end do
            !
            !==> case m = 1 through nlat-1
            !
            if (mmax < 2) then
                return
            end if

            do mp1=2, mmax
                m = mp1-1
                mb = m*(nlat-1)-(m*(m-1))/2
                mp2 = mp1+1

                if (mp1 > ndo1) then
                    exit
                end if

                do k=1, nt
                    do i=1, imm1
                        do np1=mp1, ndo1, 2
                            cr(mp1, np1, k) = cr(mp1, np1, k)-zv(np1+mb, i)*wo(i, 2*mp1-2, k) &
                                +zw(np1+mb, i)*ve(i, 2*mp1-1, k)
                            ci(mp1, np1, k) = ci(mp1, np1, k)-zv(np1+mb, i)*wo(i, 2*mp1-1, k) &
                                -zw(np1+mb, i)*ve(i, 2*mp1-2, k)
                        end do
                    end do
                end do

                if (mlat == 0) then
                    exit
                end if

                do k=1, nt
                    do np1=mp1, ndo1, 2
                        cr(mp1, np1, k) = cr(mp1, np1, k)+zw(np1+mb, imid)*ve(imid, 2*mp1-1, k)
                        ci(mp1, np1, k) = ci(mp1, np1, k)-zw(np1+mb, imid)*ve(imid, 2*mp1-2, k)
                    end do
                end do
            end do
            return
        case (6)
            !
            !==> case ityp=6 ,  v odd , w even
            !
            !    case m=0
            !
            do k=1, nt
                do i=1, imid
                    do np1=2, ndo2, 2
                        cr(1, np1, k) = cr(1, np1, k)-zv(np1, i)*we(i, 1, k)
                    end do
                end do
            end do

            do k=1, nt
                do i=1, imm1
                    do np1=3, ndo1, 2
                        br(1, np1, k) = br(1, np1, k)+zv(np1, i)*vo(i, 1, k)
                    end do
                end do
            end do
            !
            !==> case m = 1 through nlat-1
            !
            if (mmax < 2) then
                return
            end if

            do mp1=2, mmax
                m = mp1-1
                mb = m*(nlat-1)-(m*(m-1))/2
                mp2 = mp1+1

                if (mp1 > ndo1) then
                    go to 617
                end if

                do  k=1, nt
                    do  i=1, imm1
                        do  np1=mp1, ndo1, 2
                            br(mp1, np1, k) = br(mp1, np1, k)+zv(np1+mb, i)*vo(i, 2*mp1-2, k) &
                                +zw(np1+mb, i)*we(i, 2*mp1-1, k)
                            bi(mp1, np1, k) = bi(mp1, np1, k)+zv(np1+mb, i)*vo(i, 2*mp1-1, k) &
                                -zw(np1+mb, i)*we(i, 2*mp1-2, k)
                        end do
                    end do
                end do

                if (mlat == 0) then
                    go to 617
                end if

                do  k=1, nt
                    do  np1=mp1, ndo1, 2
                        br(mp1, np1, k) = br(mp1, np1, k)+zw(np1+mb, imid)*we(imid, 2*mp1-1, k)
                        bi(mp1, np1, k) = bi(mp1, np1, k)-zw(np1+mb, imid)*we(imid, 2*mp1-2, k)
                    end do
                end do

617             if (mp2 > ndo2) then
                    exit
                end if

                do k=1, nt
                    do  i=1, imm1
                        do  np1=mp2, ndo2, 2
                            cr(mp1, np1, k) = cr(mp1, np1, k)-zv(np1+mb, i)*we(i, 2*mp1-2, k) &
                                +zw(np1+mb, i)*vo(i, 2*mp1-1, k)
                            ci(mp1, np1, k) = ci(mp1, np1, k)-zv(np1+mb, i)*we(i, 2*mp1-1, k) &
                                -zw(np1+mb, i)*vo(i, 2*mp1-2, k)
                        end do
                    end do
                end do

                if (mlat == 0) then
                    exit
                end if

                do k=1, nt
                    do np1=mp2, ndo2, 2
                        cr(mp1, np1, k) = cr(mp1, np1, k)-zv(np1+mb, imid)*we(imid, 2*mp1-2, k)
                        ci(mp1, np1, k) = ci(mp1, np1, k)-zv(np1+mb, imid)*we(imid, 2*mp1-1, k)
                    end do
                end do
            end do
            return
        case (7)
            !
            !==> case ityp=7   v odd, w even, and cr and ci equal zero
            !
            !    case m=0
            !
            do k=1, nt
                do i=1, imm1
                    do np1=3, ndo1, 2
                        br(1, np1, k) = br(1, np1, k)+zv(np1, i)*vo(i, 1, k)
                    end do
                end do
            end do
            !
            !==> case m = 1 through nlat-1
            !
            if (mmax < 2) then
                return
            end if

            do mp1=2, mmax
                m = mp1-1
                mb = m*(nlat-1)-(m*(m-1))/2
                mp2 = mp1+1

                if (mp1 > ndo1) then
                    exit
                end if

                do  k=1, nt
                    do  i=1, imm1
                        do  np1=mp1, ndo1, 2
                            br(mp1, np1, k) = br(mp1, np1, k)+zv(np1+mb, i)*vo(i, 2*mp1-2, k) &
                                +zw(np1+mb, i)*we(i, 2*mp1-1, k)
                            bi(mp1, np1, k) = bi(mp1, np1, k)+zv(np1+mb, i)*vo(i, 2*mp1-1, k) &
                                -zw(np1+mb, i)*we(i, 2*mp1-2, k)
                        end do
                    end do
                end do

                if (mlat == 0) then
                    exit
                end if

                do k=1, nt
                    do np1=mp1, ndo1, 2
                        br(mp1, np1, k) = br(mp1, np1, k)+zw(np1+mb, imid)*we(imid, 2*mp1-1, k)
                        bi(mp1, np1, k) = bi(mp1, np1, k)-zw(np1+mb, imid)*we(imid, 2*mp1-2, k)
                    end do
                end do
            end do
            return
        case (8)
            !
            !==> case ityp=8   v odd, w even, and both br and bi equal zero
            !
            !    case m=0
            !
            do k=1, nt
                do i=1, imid
                    do np1=2, ndo2, 2
                        cr(1, np1, k) = cr(1, np1, k)-zv(np1, i)*we(i, 1, k)
                    end do
                end do
            end do
            !
            !     case m = 1 through nlat-1
            !
            if (mmax < 2) then
                return
            end if

            do mp1=2, mmax
                m = mp1-1
                mb = m*(nlat-1)-(m*(m-1))/2
                mp2 = mp1+1

                if (mp2 > ndo2) then
                    exit
                end if

                do k=1, nt
                    do i=1, imm1
                        do np1=mp2, ndo2, 2
                            cr(mp1, np1, k) = cr(mp1, np1, k)-zv(np1+mb, i)*we(i, 2*mp1-2, k) &
                                +zw(np1+mb, i)*vo(i, 2*mp1-1, k)
                            ci(mp1, np1, k) = ci(mp1, np1, k)-zv(np1+mb, i)*we(i, 2*mp1-1, k) &
                                -zw(np1+mb, i)*vo(i, 2*mp1-2, k)
                        end do
                    end do
                end do

                if (mlat == 0) then
                    exit
                end if

                do k=1, nt
                    do np1=mp2, ndo2, 2
                        cr(mp1, np1, k) = cr(mp1, np1, k)-zv(np1+mb, imid)*we(imid, 2*mp1-2, k)
                        ci(mp1, np1, k) = ci(mp1, np1, k)-zv(np1+mb, imid)*we(imid, 2*mp1-1, k)
                    end do
                end do
            end do
    end select

end subroutine vhaes1



subroutine vhaesi(nlat, nlon, wvhaes, lvhaes, work, lwork, dwork, &
    ldwork, ierror)
    !
    !     dwork must be of length at least 2*(nlat+1)
    !
    dimension wvhaes(lvhaes), work(lwork)
    real dwork(ldwork)
    ierror = 1
    if (nlat < 3) return
    ierror = 2
    if (nlon < 1) return
    ierror = 3
    mmax = min(nlat, (nlon+1)/2)
    imid = (nlat+1)/2
    lzimn = (imid*mmax*(nlat+nlat-mmax+1))/2
    if (lvhaes < lzimn+lzimn+nlon+15) return
    ierror = 4
    labc = 3*(max(mmax-2, 0)*(nlat+nlat-mmax-1))/2
    if (lwork < 5*nlat*imid+labc) return
    ierror = 5
    if (ldwork < 2*(nlat+1)) return
    ierror = 0
    iw1 = 3*nlat*imid+1
    idz = (mmax*(nlat+nlat-mmax+1))/2

    call vea1(nlat, nlon, imid, wvhaes, wvhaes(lzimn+1), idz, &
        work, work(iw1), dwork)

    call hrffti(nlon, wvhaes(2*lzimn+1))

end subroutine vhaesi



subroutine vea1(nlat, nlon, imid, zv, zw, idz, zin, wzvin, dwork)
    dimension zv(idz, 1), zw(idz, 1), zin(imid, nlat, 3), wzvin(1)
    real dwork(*)

    associate( mmax => min(nlat, (nlon+1)/2) )

        call zvinit(nlat, nlon, wzvin, dwork)

        do mp1=1, mmax
            m = mp1-1
            call zvin(0, nlat, nlon, m, zin, i3, wzvin)
            do np1=mp1, nlat
                mn = m*(nlat-1)-(m*(m-1))/2+np1
                zv(mn,1:imid) = zin(1:imid,np1,i3)
            end do
        end do

        call zwinit(nlat, nlon, wzvin, dwork)

        do mp1=1, mmax
            m = mp1-1
            call zwin(0, nlat, nlon, m, zin, i3, wzvin)
            do np1=mp1, nlat
                mn = m*(nlat-1)-(m*(m-1))/2+np1
                zw(mn,1:imid) = zin(1:imid,np1,i3)
            end do
        end do
    end associate

end subroutine vea1
