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
! ... file vtses.f
!
!     this file includes documentation and code for
!     subroutines vtses and vtsesi
!
! ... files which must be loaded with vtses.f
!
!     sphcom.f, hrfft.f, vhaes.f, vhses.f
!   
!   
!     subroutine vtses(nlat, nlon, ityp, nt, vt, wt, idvw, jdvw, br, bi, cr, ci, 
!    +                 mdab, ndab, wvts, lwvts, work, lwork, ierror)
!
!     given the vector harmonic analysis br, bi, cr, and ci (computed
!     by subroutine vhaes) of some vector function (v, w), this 
!     subroutine computes the vector function (vt, wt) which is
!     the derivative of (v, w) with respect to colatitude theta. vtses
!     is similar to vhses except the vector harmonics are replaced by 
!     their derivative with respect to colatitude with the result that
!     (vt, wt) is computed instead of (v, w). vt(i, j) is the derivative
!     of the colatitudinal component v(i, j) at the point theta(i) =
!     (i-1)*pi/(nlat-1) and longitude phi(j) = (j-1)*2*pi/nlon. the
!     spectral representation of (vt, wt) is given below at output
!     parameters vt, wt.
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
!                 is performed on the entire sphere. i.e. the arrays
!                 vt(i, j), wt(i, j) are computed for i=1, ..., nlat and
!                 j=1, ..., nlon.   
!
!            = 1  no symmetries exist about the equator however the
!                 the coefficients cr and ci are zero. the synthesis
!                 is performed on the entire sphere.  i.e. the arrays
!                 vt(i, j), wt(i, j) are computed for i=1, ..., nlat and
!                 j=1, ..., nlon.
!
!            = 2  no symmetries exist about the equator however the
!                 the coefficients br and bi are zero. the synthesis
!                 is performed on the entire sphere.  i.e. the arrays
!                 vt(i, j), wt(i, j) are computed for i=1, ..., nlat and 
!                 j=1, ..., nlon.
!
!            = 3  vt is odd and wt is even about the equator. the 
!                 synthesis is performed on the northern hemisphere
!                 only.  i.e., if nlat is odd the arrays vt(i, j), wt(i, j)
!                 are computed for i=1, ...,(nlat+1)/2 and j=1, ..., nlon.
!                 if nlat is even the arrays vt(i, j), wt(i, j) are computed 
!                 for i=1, ..., nlat/2 and j=1, ..., nlon.
!
!            = 4  vt is odd and wt is even about the equator and the 
!                 coefficients cr and ci are zero. the synthesis is
!                 performed on the northern hemisphere only. i.e. if
!                 nlat is odd the arrays vt(i, j), wt(i, j) are computed
!                 for i=1, ...,(nlat+1)/2 and j=1, ..., nlon. if nlat is
!                 even the arrays vt(i, j), wt(i, j) are computed for 
!                 i=1, ..., nlat/2 and j=1, ..., nlon.
!
!            = 5  vt is odd and wt is even about the equator and the 
!                 coefficients br and bi are zero. the synthesis is
!                 performed on the northern hemisphere only. i.e. if
!                 nlat is odd the arrays vt(i, j), wt(i, j) are computed
!                 for i=1, ...,(nlat+1)/2 and j=1, ..., nlon. if nlat is
!                 even the arrays vt(i, j), wt(i, j) are computed for 
!                 i=1, ..., nlat/2 and j=1, ..., nlon.
!
!            = 6  vt is even and wt is odd about the equator. the 
!                 synthesis is performed on the northern hemisphere
!                 only.  i.e., if nlat is odd the arrays vt(i, j), wt(i, j)
!                 are computed for i=1, ...,(nlat+1)/2 and j=1, ..., nlon.
!                 if nlat is even the arrays vt(i, j), wt(i, j) are computed 
!                 for i=1, ..., nlat/2 and j=1, ..., nlon.
!
!            = 7  vt is even and wt is odd about the equator and the 
!                 coefficients cr and ci are zero. the synthesis is
!                 performed on the northern hemisphere only. i.e. if
!                 nlat is odd the arrays vt(i, j), wt(i, j) are computed
!                 for i=1, ...,(nlat+1)/2 and j=1, ..., nlon. if nlat is
!                 even the arrays vt(i, j), wt(i, j) are computed for 
!                 i=1, ..., nlat/2 and j=1, ..., nlon.
!
!            = 8  vt is even and wt is odd about the equator and the 
!                 coefficients br and bi are zero. the synthesis is
!                 performed on the northern hemisphere only. i.e. if
!                 nlat is odd the arrays vt(i, j), wt(i, j) are computed
!                 for i=1, ...,(nlat+1)/2 and j=1, ..., nlon. if nlat is
!                 even the arrays vt(i, j), wt(i, j) are computed for 
!                 i=1, ..., nlat/2 and j=1, ..., nlon.
!
!     nt     the number of syntheses.  in the program that calls vtses, 
!            the arrays vt, wt, br, bi, cr, and ci can be three dimensional
!            in which case multiple syntheses will be performed.
!            the third index is the synthesis index which assumes the 
!            values k=1, ..., nt.  for a single synthesis set nt=1. the
!            discription of the remaining parameters is simplified
!            by assuming that nt=1 or that all the arrays are two
!            dimensional.
!
!     idvw   the first dimension of the arrays vt, wt as it appears in
!            the program that calls vtses. if ityp .le. 2 then idvw
!            must be at least nlat.  if ityp .gt. 2 and nlat is
!            even then idvw must be at least nlat/2. if ityp .gt. 2
!            and nlat is odd then idvw must be at least(nlat+1)/2.
!
!     jdvw   the second dimension of the arrays vt, wt as it appears in
!            the program that calls vtses. jdvw must be at least nlon.
!
!     br, bi  two or three dimensional arrays (see input parameter nt)
!     cr, ci  that contain the vector spherical harmonic coefficients
!            of (v, w) as computed by subroutine vhaes.
!            
!     mdab   the first dimension of the arrays br, bi, cr, and ci as it
!            appears in the program that calls vtses. mdab must be at
!            least min(nlat, nlon/2) if nlon is even or at least
!            min(nlat, (nlon+1)/2) if nlon is odd.
!
!     ndab   the second dimension of the arrays br, bi, cr, and ci as it
!            appears in the program that calls vtses. ndab must be at
!            least nlat.
!
!     wvts   an array which must be initialized by subroutine vtsesi.
!            once initialized, wvts can be used repeatedly by vtses
!            as long as nlon and nlat remain unchanged.  wvts must
!            not be altered between calls of vtses.
!
!     lwvts  the dimension of the array wvts as it appears in the
!            program that calls vtses. define
!
!               l1 = min(nlat, nlon/2) if nlon is even or
!               l1 = min(nlat, (nlon+1)/2) if nlon is odd
!
!            and
!
!               l2 = nlat/2        if nlat is even or
!               l2 =(nlat+1)/2    if nlat is odd
!
!            then lwvts must be at least
!
!                 l1*l2*(nlat+nlat-l1+1)+nlon+15
!
!
!     work   a work array that does not have to be saved.
!
!     lwork  the dimension of the array work as it appears in the
!            program that calls vtses. define
!
!               l2 = nlat/2        if nlat is even or
!               l2 =(nlat+1)/2    if nlat is odd
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
!     vt, wt  two or three dimensional arrays (see input parameter nt)
!            in which the derivative of (v, w) with respect to 
!            colatitude theta is stored. vt(i, j), wt(i, j) contain the
!            derivatives at colatitude theta(i) = (i-1)*pi/(nlat-1)
!            and longitude phi(j) = (j-1)*2*pi/nlon. the index ranges
!            are defined above at the input parameter ityp. vt and wt
!            are computed from the formulas for v and w given in 
!            subroutine vhses but with vbar and wbar replaced with
!            their derivatives with respect to colatitude. these
!            derivatives are denoted by vtbar and wtbar. 
!
!
!   *************************************************************
!
!   in terms of real variables this expansion takes the form
!
!             for i=1, ..., nlat and  j=1, ..., nlon
!
!     vt(i, j) = the sum from n=1 to n=nlat-1 of
!
!               .5*br(1, n+1)*vtbar(0, n, theta(i))
!
!     plus the sum from m=1 to m=mmax-1 of the sum from n=m to 
!     n=nlat-1 of the real part of
!
!       (br(m+1, n+1)*vtbar(m, n, theta(i))
!                   -ci(m+1, n+1)*wtbar(m, n, theta(i)))*cos(m*phi(j))
!      -(bi(m+1, n+1)*vtbar(m, n, theta(i))
!                   +cr(m+1, n+1)*wtbar(m, n, theta(i)))*sin(m*phi(j))
!
!    and for i=1, ..., nlat and  j=1, ..., nlon
!
!     wt(i, j) = the sum from n=1 to n=nlat-1 of
!
!              -.5*cr(1, n+1)*vtbar(0, n, theta(i))
!
!     plus the sum from m=1 to m=mmax-1 of the sum from n=m to
!     n=nlat-1 of the real part of
!
!      -(cr(m+1, n+1)*vtbar(m, n, theta(i))
!                   +bi(m+1, n+1)*wtbar(m, n, theta(i)))*cos(m*phi(j))
!      +(ci(m+1, n+1)*vtbar(m, n, theta(i))
!                   -br(m+1, n+1)*wtbar(m, n, theta(i)))*sin(m*phi(j))
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
!            = 9  error in the specification of lwvts
!            = 10 error in the specification of lwork
!
!
! *******************************************************************
!
!     subroutine vtsesi(nlat, nlon, wvts, lwvts, work, lwork, dwork, ldwork, 
!    +                  ierror)
!
!     subroutine vtsesi initializes the array wvts which can then be
!     used repeatedly by subroutine vtses until nlat or nlon is changed.
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
!     lwvts  the dimension of the array wvts as it appears in the
!            program that calls vtses. define
!
!               l1 = min(nlat, nlon/2) if nlon is even or
!               l1 = min(nlat, (nlon+1)/2) if nlon is odd
!
!            and
!
!               l2 = nlat/2        if nlat is even or
!               l2 =(nlat+1)/2    if nlat is odd
!
!            then lwvts must be at least
!
!                  l1*l2*(nlat+nlat-l1+1)+nlon+15
!
!
!     work   a work array that does not have to be saved.
!
!     lwork  the dimension of the array work as it appears in the
!            program that calls vtses. lwork must be at least
!
!              3*(max(l1-2, 0)*(nlat+nlat-l1-1))/2+5*l2*nlat
!
!     dwork  a real work array that does have to be saved.
!
!     ldwork the length of dwork.  ldwork must be at least 2*(nlat+1)
!
!     **************************************************************
!
!     output parameters
!
!     wvts   an array which is initialized for use by subroutine vtses.
!            once initialized, wvts can be used repeatedly by vtses
!            as long as nlat or nlon remain unchanged.  wvts must not
!            be altered between calls of vtses.
!
!
!     ierror = 0  no errors
!            = 1  error in the specification of nlat
!            = 2  error in the specification of nlon
!            = 3  error in the specification of lwvts
!            = 4  error in the specification of lwork
!            = 5  error in the specification of ldwork
!
!
subroutine vtses(nlat, nlon, ityp, nt, vt, wt, idvw, jdvw, br, bi, cr, ci, &
    mdab, ndab, wvts, lwvts, work, lwork, ierror)
    !
    dimension vt(idvw, jdvw, 1), wt(idvw, jdvw, 1), br(mdab, ndab, 1), &
        bi(mdab, ndab, 1), cr(mdab, ndab, 1), ci(mdab, ndab, 1), &
        work(1), wvts(1)

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
    idz = (mmax*(nlat+nlat-mmax+1))/2
    lzimn = idz*imid
    if (lwvts < lzimn+lzimn+nlon+15) return
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

    call vtses1(nlat, nlon, ityp, nt, imid, idvw, jdvw, vt, wt, mdab, ndab, &
        br, bi, cr, ci, idv, work, work(iw1), work(iw2), work(iw3), &
        work(iw4), idz, wvts, wvts(jw1), wvts(jw2))

end subroutine vtses



subroutine vtses1(nlat, nlon, ityp, nt, imid, idvw, jdvw, vt, wt, mdab, &
    ndab, br, bi, cr, ci, idv, vte, vto, wte, wto, work, idz, vb, wb, wrfft)
    dimension vt(idvw, jdvw, 1), wt(idvw, jdvw, 1), br(mdab, ndab, 1), &
        bi(mdab, ndab, 1), cr(mdab, ndab, 1), ci(mdab, ndab, 1), &
        vte(idv, nlon, 1), vto(idv, nlon, 1), wte(idv, nlon, 1), &
        wto(idv, nlon, 1), work(1), wrfft(1), &
        vb(imid, 1), wb(imid, 1)

    integer :: i, j, k, hack_counter !! Counters

    nlp1 = nlat+1
    mlat = mod(nlat, 2)
    mlon = mod(nlon, 2)
    mmax = min(nlat, (nlon+1)/2)

    if (mlat /= 0) then
        imm1 = imid-1
        ndo1 = nlat-1
        ndo2 = nlat
    else
        imm1 = imid
        ndo1 = nlat
        ndo2 = nlat-1
    end if
    !
    !==> Initialize
    !
    vte = 0.0
    wte = 0.0

    !
    !==> Remark:
    !
    !    Fortran 2008:
    !    The EXIT statement provides one way of terminating
    !    a loop, or completing execution of another construct.
    !    Unfortunately, many compilers don't support this feature yet.
    !
    !    GCC's gfortran 5.1 doesn't support EXIT within the SELECT CASE
    !    construct yet.
    !
    !    We exit a DO loop labeled "dumb_hack" to avoid the dreaded GO TO
    !
    dumb_hack: do hack_counter=1,1
        select case (ityp)
            case (0)
                !
                !     case ityp=0   no symmetries
                !
                !     case m = 0
                !
                do k=1, nt
                    do np1=2, ndo2, 2
                        do i=1, imm1
                            vto(i, 1, k)=vto(i, 1, k)+br(1, np1, k)*vb(i, np1)
                            wto(i, 1, k)=wto(i, 1, k)-cr(1, np1, k)*vb(i, np1)
                        end do
                    end do
                end do
                do k=1, nt
                    do np1=3, ndo1, 2
                        do i=1, imid
                            vte(i, 1, k)=vte(i, 1, k)+br(1, np1, k)*vb(i, np1)
                            wte(i, 1, k)=wte(i, 1, k)-cr(1, np1, k)*vb(i, np1)
                        end do
                    end do
                end do
                !
                !     case m = 1 through nlat-1
                !
                if (mmax < 2) then
                    exit dumb_hack
                end if

                do mp1=2, mmax
                    m = mp1-1
                    mb = m*(nlat-1)-(m*(m-1))/2
                    mp2 = mp1+1
                    if (mp1 > ndo1) then
                        if (mp2 > ndo2) then
                            exit
                        end if
                    end if
                    do k=1, nt
                        do np1=mp1, ndo1, 2
                            mn = mb+np1
                            do i=1, imm1
                                vte(i, 2*mp1-2, k) = vte(i, 2*mp1-2, k)+br(mp1, np1, k)*vb(i, mn)
                                vto(i, 2*mp1-2, k) = vto(i, 2*mp1-2, k)-ci(mp1, np1, k)*wb(i, mn)
                                vte(i, 2*mp1-1, k) = vte(i, 2*mp1-1, k)+bi(mp1, np1, k)*vb(i, mn)
                                vto(i, 2*mp1-1, k) = vto(i, 2*mp1-1, k)+cr(mp1, np1, k)*wb(i, mn)
                                wte(i, 2*mp1-2, k) = wte(i, 2*mp1-2, k)-cr(mp1, np1, k)*vb(i, mn)
                                wto(i, 2*mp1-2, k) = wto(i, 2*mp1-2, k)-bi(mp1, np1, k)*wb(i, mn)
                                wte(i, 2*mp1-1, k) = wte(i, 2*mp1-1, k)-ci(mp1, np1, k)*vb(i, mn)
                                wto(i, 2*mp1-1, k) = wto(i, 2*mp1-1, k)+br(mp1, np1, k)*wb(i, mn)
                            end do
                            if (mlat == 0) then
                                exit
                            end if

                            vte(imid, 2*mp1-2, k) = vte(imid, 2*mp1-2, k) &
                                +br(mp1, np1, k)*vb(imid, mn)
                            vte(imid, 2*mp1-1, k) = vte(imid, 2*mp1-1, k) &
                                +bi(mp1, np1, k)*vb(imid, mn)
                            wte(imid, 2*mp1-2, k) = wte(imid, 2*mp1-2, k) &
                                -cr(mp1, np1, k)*vb(imid, mn)
                            wte(imid, 2*mp1-1, k) = wte(imid, 2*mp1-1, k) &
                                -ci(mp1, np1, k)*vb(imid, mn)
                        end do
                    end do

                    do k=1, nt
                        do np1=mp2, ndo2, 2
                            mn = mb+np1
                            do i=1, imm1
                                vto(i, 2*mp1-2, k) = vto(i, 2*mp1-2, k)+br(mp1, np1, k)*vb(i, mn)
                                vte(i, 2*mp1-2, k) = vte(i, 2*mp1-2, k)-ci(mp1, np1, k)*wb(i, mn)
                                vto(i, 2*mp1-1, k) = vto(i, 2*mp1-1, k)+bi(mp1, np1, k)*vb(i, mn)
                                vte(i, 2*mp1-1, k) = vte(i, 2*mp1-1, k)+cr(mp1, np1, k)*wb(i, mn)
                                wto(i, 2*mp1-2, k) = wto(i, 2*mp1-2, k)-cr(mp1, np1, k)*vb(i, mn)
                                wte(i, 2*mp1-2, k) = wte(i, 2*mp1-2, k)-bi(mp1, np1, k)*wb(i, mn)
                                wto(i, 2*mp1-1, k) = wto(i, 2*mp1-1, k)-ci(mp1, np1, k)*vb(i, mn)
                                wte(i, 2*mp1-1, k) = wte(i, 2*mp1-1, k)+br(mp1, np1, k)*wb(i, mn)
                            end do
                            if (mlat == 0) then
                                exit
                            end if
                            vte(imid, 2*mp1-2, k) = vte(imid, 2*mp1-2, k) &
                                -ci(mp1, np1, k)*wb(imid, mn)
                            vte(imid, 2*mp1-1, k) = vte(imid, 2*mp1-1, k) &
                                +cr(mp1, np1, k)*wb(imid, mn)
                            wte(imid, 2*mp1-2, k) = wte(imid, 2*mp1-2, k) &
                                -bi(mp1, np1, k)*wb(imid, mn)
                            wte(imid, 2*mp1-1, k) = wte(imid, 2*mp1-1, k) &
                                +br(mp1, np1, k)*wb(imid, mn)
                        end do
                    end do
                end do
            case (1)
                !
                !==> case ityp=1   no symmetries,  cr and ci equal zero
                !
                !     case m = 0
                !
                do k=1, nt
                    do np1=2, ndo2, 2
                        do i=1, imm1
                            vto(i, 1, k)=vto(i, 1, k)+br(1, np1, k)*vb(i, np1)
                        end do
                    end do
                end do

                do k=1, nt
                    do np1=3, ndo1, 2
                        do i=1, imid
                            vte(i, 1, k)=vte(i, 1, k)+br(1, np1, k)*vb(i, np1)
                        end do
                    end do
                end do
                !
                !     case m = 1 through nlat-1
                !
                if (mmax < 2) then
                    exit dumb_hack
                end if

                do mp1=2, mmax
                    m = mp1-1
                    mb = m*(nlat-1)-(m*(m-1))/2
                    mp2 = mp1+1

                    if (mp1 > ndo1) then
                        if (mp2 > ndo2) then
                            exit
                        end if
                    end if

                    do k=1, nt
                        do np1=mp1, ndo1, 2
                            mn = mb+np1
                            do i=1, imm1
                                vte(i, 2*mp1-2, k) = vte(i, 2*mp1-2, k)+br(mp1, np1, k)*vb(i, mn)
                                vte(i, 2*mp1-1, k) = vte(i, 2*mp1-1, k)+bi(mp1, np1, k)*vb(i, mn)
                                wto(i, 2*mp1-2, k) = wto(i, 2*mp1-2, k)-bi(mp1, np1, k)*wb(i, mn)
                                wto(i, 2*mp1-1, k) = wto(i, 2*mp1-1, k)+br(mp1, np1, k)*wb(i, mn)
                            end do
                            if (mlat == 0) then
                            end if

                            vte(imid, 2*mp1-2, k) = vte(imid, 2*mp1-2, k) &
                                +br(mp1, np1, k)*vb(imid, mn)
                            vte(imid, 2*mp1-1, k) = vte(imid, 2*mp1-1, k) &
                                +bi(mp1, np1, k)*vb(imid, mn)
                        end do
                    end do


                    do k=1, nt
                        do np1=mp2, ndo2, 2
                            mn = mb+np1
                            do i=1, imm1
                                vto(i, 2*mp1-2, k) = vto(i, 2*mp1-2, k)+br(mp1, np1, k)*vb(i, mn)
                                vto(i, 2*mp1-1, k) = vto(i, 2*mp1-1, k)+bi(mp1, np1, k)*vb(i, mn)
                                wte(i, 2*mp1-2, k) = wte(i, 2*mp1-2, k)-bi(mp1, np1, k)*wb(i, mn)
                                wte(i, 2*mp1-1, k) = wte(i, 2*mp1-1, k)+br(mp1, np1, k)*wb(i, mn)
                            end do
                            if (mlat == 0) then
                                exit
                            end if
                            wte(imid, 2*mp1-2, k) = wte(imid, 2*mp1-2, k) &
                                -bi(mp1, np1, k)*wb(imid, mn)
                            wte(imid, 2*mp1-1, k) = wte(imid, 2*mp1-1, k) &
                                +br(mp1, np1, k)*wb(imid, mn)
                        end do
                    end do
                end do
            case (2)
                !
                !     case ityp=2   no symmetries,  br and bi are equal to zero
                !
                !     case m = 0
                !
                do k=1, nt
                    do np1=2, ndo2, 2
                        do i=1, imm1
                            wto(i, 1, k)=wto(i, 1, k)-cr(1, np1, k)*vb(i, np1)
                        end do
                    end do
                end do

                do k=1, nt
                    do np1=3, ndo1, 2
                        do i=1, imid
                            wte(i, 1, k)=wte(i, 1, k)-cr(1, np1, k)*vb(i, np1)
                        end do
                    end do
                end do
                !
                !     case m = 1 through nlat-1
                !
                if (mmax < 2) then
                    exit dumb_hack
                end if

                do mp1=2, mmax
                    m = mp1-1
                    mb = m*(nlat-1)-(m*(m-1))/2
                    mp2 = mp1+1
                    if (mp1 > ndo1) then
                        if (mp2 > ndo2) then
                            exit
                        end if
                    end if

                    do k=1, nt
                        do np1=mp1, ndo1, 2
                            mn = mb+np1
                            do i=1, imm1
                                vto(i, 2*mp1-2, k) = vto(i, 2*mp1-2, k)-ci(mp1, np1, k)*wb(i, mn)
                                vto(i, 2*mp1-1, k) = vto(i, 2*mp1-1, k)+cr(mp1, np1, k)*wb(i, mn)
                                wte(i, 2*mp1-2, k) = wte(i, 2*mp1-2, k)-cr(mp1, np1, k)*vb(i, mn)
                                wte(i, 2*mp1-1, k) = wte(i, 2*mp1-1, k)-ci(mp1, np1, k)*vb(i, mn)
                            end do
                            if (mlat == 0) then
                                exit
                            end if
                            wte(imid, 2*mp1-2, k) = wte(imid, 2*mp1-2, k) &
                                -cr(mp1, np1, k)*vb(imid, mn)
                            wte(imid, 2*mp1-1, k) = wte(imid, 2*mp1-1, k) &
                                -ci(mp1, np1, k)*vb(imid, mn)
                        end do
                    end do

                    do k=1, nt
                        do np1=mp2, ndo2, 2
                            mn = mb+np1
                            do i=1, imm1
                                vte(i, 2*mp1-2, k) = vte(i, 2*mp1-2, k)-ci(mp1, np1, k)*wb(i, mn)
                                vte(i, 2*mp1-1, k) = vte(i, 2*mp1-1, k)+cr(mp1, np1, k)*wb(i, mn)
                                wto(i, 2*mp1-2, k) = wto(i, 2*mp1-2, k)-cr(mp1, np1, k)*vb(i, mn)
                                wto(i, 2*mp1-1, k) = wto(i, 2*mp1-1, k)-ci(mp1, np1, k)*vb(i, mn)
                            end do
                            if (mlat == 0) then
                                exit
                            end if
                            vte(imid, 2*mp1-2, k) = vte(imid, 2*mp1-2, k) &
                                -ci(mp1, np1, k)*wb(imid, mn)
                            vte(imid, 2*mp1-1, k) = vte(imid, 2*mp1-1, k) &
                                +cr(mp1, np1, k)*wb(imid, mn)
                        end do
                    end do
                end do
            case (3)
                 !
                 !==> case ityp=3   v odd,  w even
                 !
                 !     case m = 0
                 !
                do k=1, nt
                    do np1=2, ndo2, 2
                        do i=1, imm1
                            vto(i, 1, k)=vto(i, 1, k)+br(1, np1, k)*vb(i, np1)
                        end do
                    end do
                end do

                do k=1, nt
                    do np1=3, ndo1, 2
                        do i=1, imid
                            wte(i, 1, k)=wte(i, 1, k)-cr(1, np1, k)*vb(i, np1)
                        end do
                    end do
                end do
                !
                !     case m = 1 through nlat-1
                !
                if (mmax < 2) then
                    exit dumb_hack
                end if

                do mp1=2, mmax
                    m = mp1-1
                    mb = m*(nlat-1)-(m*(m-1))/2
                    mp2 = mp1+1
                    if (mp1 > ndo1) then
                        if (mp2 > ndo2) then
                            exit
                        end if

                    end if

                    do k=1, nt
                        do np1=mp1, ndo1, 2
                            mn = mb+np1
                            do i=1, imm1
                                vto(i, 2*mp1-2, k) = vto(i, 2*mp1-2, k)-ci(mp1, np1, k)*wb(i, mn)
                                vto(i, 2*mp1-1, k) = vto(i, 2*mp1-1, k)+cr(mp1, np1, k)*wb(i, mn)
                                wte(i, 2*mp1-2, k) = wte(i, 2*mp1-2, k)-cr(mp1, np1, k)*vb(i, mn)
                                wte(i, 2*mp1-1, k) = wte(i, 2*mp1-1, k)-ci(mp1, np1, k)*vb(i, mn)
                            end do

                            if (mlat == 0) then
                                exit
                            end if

                            wte(imid, 2*mp1-2, k) = wte(imid, 2*mp1-2, k) &
                                -cr(mp1, np1, k)*vb(imid, mn)
                            wte(imid, 2*mp1-1, k) = wte(imid, 2*mp1-1, k) &
                                -ci(mp1, np1, k)*vb(imid, mn)
                        end do
                    end do

                    do k=1, nt
                        do np1=mp2, ndo2, 2
                            mn = mb+np1
                            do i=1, imm1
                                vto(i, 2*mp1-2, k) = vto(i, 2*mp1-2, k)+br(mp1, np1, k)*vb(i, mn)
                                vto(i, 2*mp1-1, k) = vto(i, 2*mp1-1, k)+bi(mp1, np1, k)*vb(i, mn)
                                wte(i, 2*mp1-2, k) = wte(i, 2*mp1-2, k)-bi(mp1, np1, k)*wb(i, mn)
                                wte(i, 2*mp1-1, k) = wte(i, 2*mp1-1, k)+br(mp1, np1, k)*wb(i, mn)
                            end do
                            if (mlat == 0) then
                                exit
                            end if
                            wte(imid, 2*mp1-2, k) = wte(imid, 2*mp1-2, k) &
                                -bi(mp1, np1, k)*wb(imid, mn)
                            wte(imid, 2*mp1-1, k) = wte(imid, 2*mp1-1, k) &
                                +br(mp1, np1, k)*wb(imid, mn)
                        end do
                    end do
                end do
            case (4)
                !
                !==> case ityp=4   v odd,  w even, and both cr and ci equal zero
                !
                !     case m = 0
                !
                do k=1, nt
                    do np1=2, ndo2, 2
                        do i=1, imm1
                            vto(i, 1, k)=vto(i, 1, k)+br(1, np1, k)*vb(i, np1)
                        end do
                    end do
                end do
                !
                !     case m = 1 through nlat-1
                !
                if (mmax < 2) then
                    exit dumb_hack
                end if

                do mp1=2, mmax
                    m = mp1-1
                    mb = m*(nlat-1)-(m*(m-1))/2
                    mp2 = mp1+1

                    if (mp2 > ndo2) then
                        exit
                    end if

                    do k=1, nt
                        do np1=mp2, ndo2, 2
                            mn = mb+np1
                            do i=1, imm1
                                vto(i, 2*mp1-2, k) = vto(i, 2*mp1-2, k)+br(mp1, np1, k)*vb(i, mn)
                                vto(i, 2*mp1-1, k) = vto(i, 2*mp1-1, k)+bi(mp1, np1, k)*vb(i, mn)
                                wte(i, 2*mp1-2, k) = wte(i, 2*mp1-2, k)-bi(mp1, np1, k)*wb(i, mn)
                                wte(i, 2*mp1-1, k) = wte(i, 2*mp1-1, k)+br(mp1, np1, k)*wb(i, mn)
                            end do
                            if (mlat == 0) then
                                exit
                            end if
                            wte(imid, 2*mp1-2, k) = wte(imid, 2*mp1-2, k) &
                                -bi(mp1, np1, k)*wb(imid, mn)
                            wte(imid, 2*mp1-1, k) = wte(imid, 2*mp1-1, k) &
                                +br(mp1, np1, k)*wb(imid, mn)
                        end do
                    end do
                end do
            case (5)
                !
                !==> case ityp=5   v odd,  w even,     br and bi equal zero
                !
                !     case m = 0
                !
                do k=1, nt
                    do np1=3, ndo1, 2
                        do i=1, imid
                            wte(i, 1, k)=wte(i, 1, k)-cr(1, np1, k)*vb(i, np1)
                        end do
                    end do
                end do
                !
                !     case m = 1 through nlat-1
                !
                if (mmax < 2) then
                    exit dumb_hack
                end if

                do mp1=2, mmax
                    m = mp1-1
                    mb = m*(nlat-1)-(m*(m-1))/2
                    mp2 = mp1+1
                    if (mp1 > ndo1) then
                        exit
                    end if

                    do k=1, nt
                        do np1=mp1, ndo1, 2
                            mn = mb+np1
                            do i=1, imm1
                                vto(i, 2*mp1-2, k) = vto(i, 2*mp1-2, k)-ci(mp1, np1, k)*wb(i, mn)
                                vto(i, 2*mp1-1, k) = vto(i, 2*mp1-1, k)+cr(mp1, np1, k)*wb(i, mn)
                                wte(i, 2*mp1-2, k) = wte(i, 2*mp1-2, k)-cr(mp1, np1, k)*vb(i, mn)
                                wte(i, 2*mp1-1, k) = wte(i, 2*mp1-1, k)-ci(mp1, np1, k)*vb(i, mn)
                            end do
                            if (mlat == 0) then
                                exit
                            end if
                            wte(imid, 2*mp1-2, k) = wte(imid, 2*mp1-2, k) &
                                -cr(mp1, np1, k)*vb(imid, mn)
                            wte(imid, 2*mp1-1, k) = wte(imid, 2*mp1-1, k) &
                                -ci(mp1, np1, k)*vb(imid, mn)
                        end do
                    end do
                end do
            case (6)
                 !
                 !==> case ityp=6   v even  ,  w odd
                 !
                 !     case m = 0
                 !
                do k=1, nt
                    do np1=2, ndo2, 2
                        do i=1, imm1
                            wto(i, 1, k)=wto(i, 1, k)-cr(1, np1, k)*vb(i, np1)
                        end do
                    end do
                end do
                do k=1, nt
                    do np1=3, ndo1, 2
                        do i=1, imid
                            vte(i, 1, k)=vte(i, 1, k)+br(1, np1, k)*vb(i, np1)
                        end do
                    end do
                end do
                !
                !     case m = 1 through nlat-1
                !
                if (mmax < 2) then
                    exit dumb_hack
                end if

                do mp1=2, mmax
                    m = mp1-1
                    mb = m*(nlat-1)-(m*(m-1))/2
                    mp2 = mp1+1

                    if (mp1 > ndo1) then
                        if (mp2 > ndo2) then
                            exit
                        end if
                    end if

                    do k=1, nt
                        do np1=mp1, ndo1, 2
                            mn = mb+np1
                            do i=1, imm1
                                vte(i, 2*mp1-2, k) = vte(i, 2*mp1-2, k)+br(mp1, np1, k)*vb(i, mn)
                                vte(i, 2*mp1-1, k) = vte(i, 2*mp1-1, k)+bi(mp1, np1, k)*vb(i, mn)
                                wto(i, 2*mp1-2, k) = wto(i, 2*mp1-2, k)-bi(mp1, np1, k)*wb(i, mn)
                                wto(i, 2*mp1-1, k) = wto(i, 2*mp1-1, k)+br(mp1, np1, k)*wb(i, mn)
                            end do
                            if (mlat == 0) then
                                exit
                            end if
                            vte(imid, 2*mp1-2, k) = vte(imid, 2*mp1-2, k) &
                                +br(mp1, np1, k)*vb(imid, mn)
                            vte(imid, 2*mp1-1, k) = vte(imid, 2*mp1-1, k) &
                                +bi(mp1, np1, k)*vb(imid, mn)
                        end do
                    end do

                    do k=1, nt
                        do np1=mp2, ndo2, 2
                            mn = mb+np1
                            do i=1, imm1
                                vte(i, 2*mp1-2, k) = vte(i, 2*mp1-2, k)-ci(mp1, np1, k)*wb(i, mn)
                                vte(i, 2*mp1-1, k) = vte(i, 2*mp1-1, k)+cr(mp1, np1, k)*wb(i, mn)
                                wto(i, 2*mp1-2, k) = wto(i, 2*mp1-2, k)-cr(mp1, np1, k)*vb(i, mn)
                                wto(i, 2*mp1-1, k) = wto(i, 2*mp1-1, k)-ci(mp1, np1, k)*vb(i, mn)
                            end do
                            if (mlat == 0) then
                                exit
                            end if
                            vte(imid, 2*mp1-2, k) = vte(imid, 2*mp1-2, k) &
                                -ci(mp1, np1, k)*wb(imid, mn)
                            vte(imid, 2*mp1-1, k) = vte(imid, 2*mp1-1, k) &
                                +cr(mp1, np1, k)*wb(imid, mn)
                        end do
                    end do
                end do
            case (7)
                !
                !==> case ityp=7   v even, w odd   cr and ci equal zero
                !
                !     case m = 0
                !
                do k=1, nt
                    do np1=3, ndo1, 2
                        do i=1, imid
                            vte(i, 1, k)=vte(i, 1, k)+br(1, np1, k)*vb(i, np1)
                        end do
                    end do
                end do
                !
                !     case m = 1 through nlat-1
                !
                if (mmax < 2) then
                    exit dumb_hack
                end if

                do mp1=2, mmax
                    m = mp1-1
                    mb = m*(nlat-1)-(m*(m-1))/2
                    mp2 = mp1+1

                    if (mp1 > ndo1) then
                        exit
                    end if

                    do k=1, nt
                        do np1=mp1, ndo1, 2
                            mn = mb+np1
                            do i=1, imm1
                                vte(i, 2*mp1-2, k) = vte(i, 2*mp1-2, k)+br(mp1, np1, k)*vb(i, mn)
                                vte(i, 2*mp1-1, k) = vte(i, 2*mp1-1, k)+bi(mp1, np1, k)*vb(i, mn)
                                wto(i, 2*mp1-2, k) = wto(i, 2*mp1-2, k)-bi(mp1, np1, k)*wb(i, mn)
                                wto(i, 2*mp1-1, k) = wto(i, 2*mp1-1, k)+br(mp1, np1, k)*wb(i, mn)
                            end do
                            if (mlat == 0) then
                                exit
                            end if
                            vte(imid, 2*mp1-2, k) = vte(imid, 2*mp1-2, k) &
                                +br(mp1, np1, k)*vb(imid, mn)
                            vte(imid, 2*mp1-1, k) = vte(imid, 2*mp1-1, k) &
                                +bi(mp1, np1, k)*vb(imid, mn)
                        end do
                    end do
                end do
            case (8)
                !
                !==> case ityp=8   v even,  w odd,   br and bi equal zero
                !
                !     case m = 0
                !
                do k=1, nt
                    do np1=2, ndo2, 2
                        do i=1, imm1
                            wto(i, 1, k)=wto(i, 1, k)-cr(1, np1, k)*vb(i, np1)
                        end do
                    end do
                end do
                !
                !     case m = 1 through nlat-1
                !
                if (mmax < 2) then
                    exit dumb_hack
                end if

                do mp1=2, mmax
                    m = mp1-1
                    mb = m*(nlat-1)-(m*(m-1))/2
                    mp2 = mp1+1
                    if (mp2 > ndo2) then
                        exit
                    end if

                    do k=1, nt
                        do np1=mp2, ndo2, 2
                            mn = mb+np1
                            do i=1, imm1
                                vte(i, 2*mp1-2, k) = vte(i, 2*mp1-2, k)-ci(mp1, np1, k)*wb(i, mn)
                                vte(i, 2*mp1-1, k) = vte(i, 2*mp1-1, k)+cr(mp1, np1, k)*wb(i, mn)
                                wto(i, 2*mp1-2, k) = wto(i, 2*mp1-2, k)-cr(mp1, np1, k)*vb(i, mn)
                                wto(i, 2*mp1-1, k) = wto(i, 2*mp1-1, k)-ci(mp1, np1, k)*vb(i, mn)
                            end do
                            if (mlat == 0) then
                                exit
                            end if
                            vte(imid, 2*mp1-2, k) = vte(imid, 2*mp1-2, k) &
                                -ci(mp1, np1, k)*wb(imid, mn)
                            vte(imid, 2*mp1-1, k) = vte(imid, 2*mp1-1, k) &
                                +cr(mp1, np1, k)*wb(imid, mn)
                        end do
                    end do
                end do
        end select
    end do dumb_hack

    do k=1, nt
        call hrfftb(idv, nlon, vte(1, 1, k), idv, wrfft, work)
        call hrfftb(idv, nlon, wte(1, 1, k), idv, wrfft, work)
    end do

    select case (ityp)
        case (0,1,2)
            do i=1, imm1
                vt(i,:,:) = 0.5*( vte(i,:,:)+vto(i,:,:) )
                wt(i,:,:) = 0.5*( wte(i,:,:)+wto(i,:,:) )
                vt(nlp1-i,:,:) = 0.5*( vte(i,:,:)-vto(i, :,:) )
                wt(nlp1-i,:,:) = 0.5*( wte(i,:,:)-wto(i,:,:) )
            end do
        case default
            vt(1:imm1,:,:) = 0.5 * vte(1:imm1,:,:)
            wt(1:imm1,:,:) = 0.5 * wte(1:imm1,:,:)
    end select

    select case (mlat)
        case (0)
            return
        case default
            vt(imid,:,:) = 0.5 * vte(imid,:,:)
            wt(imid,:,:) = 0.5 * wte(imid,:,:)
    end select

end subroutine vtses1


subroutine vtsesi(nlat, nlon, wvts, lwvts, work, lwork, dwork, ldwork, &
    ierror)
    !
    dimension wvts(lwvts), work(lwork)
    real dwork(ldwork)
    ierror = 1
    if(nlat < 3) return
    ierror = 2
    if (nlon < 1) return
    ierror = 3
    mmax = min(nlat, (nlon+1)/2)
    imid =(nlat+1)/2
    lzimn = (imid*mmax*(nlat+nlat-mmax+1))/2
    if (lwvts < lzimn+lzimn+nlon+15) return
    ierror = 4
    labc = 3*(max(mmax-2, 0)*(nlat+nlat-mmax-1))/2
    if (lwork < 5*nlat*imid+labc) return
    ierror = 5
    if (ldwork < 2*(nlat+1)) return
    ierror = 0
    iw1 = 3*nlat*imid+1
    idz = (mmax*(nlat+nlat-mmax+1))/2

    call vet1(nlat, nlon, imid, wvts, wvts(lzimn+1), idz, work, work(iw1), dwork)

    call hrffti(nlon, wvts(2*lzimn+1))

end subroutine vtsesi



subroutine vet1(nlat, nlon, imid, vb, wb, idz, vin, wzvin, dwork)
    dimension vb(imid, *), wb(imid, *), vin(imid, nlat, 3), wzvin(*)
    real dwork(*)

    associate( mmax => min(nlat, (nlon+1)/2) )

        call vtinit(nlat, nlon, wzvin, dwork)

        do mp1=1, mmax
            m = mp1-1
            call vbin(0, nlat, nlon, m, vin, i3, wzvin)
            do np1=mp1, nlat
                mn = m*(nlat-1)-(m*(m-1))/2+np1
                vb(1:imid, mn) = vin(1:imid, np1, i3)
            end do
        end do

        call wtinit(nlat, nlon, wzvin, dwork)

        do mp1=1, mmax
            m = mp1-1
            call wbin(0, nlat, nlon, m, vin, i3, wzvin)
            do np1=mp1, nlat
                mn = m*(nlat-1)-(m*(m-1))/2+np1
                wb(1:imid, mn) = vin(1:imid, np1, i3)
            end do
        end do
    end associate

end subroutine vet1
