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
! ... file shigs.f
!
!     this file contains code and documentation for subroutine shigs
!
! ... files which must be loaded with shigs.f
!
!     sphcom.f, hrfft.f, gaqd.f
!
!     3/6/98
!
! *** shigs is functionally the same as shagsi or shsgsi.  It
!     is included in this version because legacy codes, using
!     older versions of spherepack called shigs to initialize
!     the saved work space wshigs for either shags or shsgs
!     Its arguments are identical to those of shagsi or shsgsi.
!
! ****************************************************************
!
!     subroutine shigs(nlat, nlon, wshigs, lshigs, work, lwork, dwork, ldwork, 
!    +                 ierror)
!
!     subroutine shigs initializes the array wshigs which can then
!     be used repeatedly by subroutines shags, shsgs. it precomputes
!     and stores in wshigs quantities such as gaussian weights, 
!     legendre polynomial coefficients, and fft trigonometric tables.
!
!     input parameters
!
!     nlat   the number of points in the gaussian colatitude grid on the
!            full sphere. these lie in the interval (0, pi) and are compu
!            in radians in theta(1), ..., theta(nlat) by subroutine gaqd.
!            if nlat is odd the equator will be included as the grid poi
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
!            than or equal to 4. the efficiency of the computation is
!            improved when nlon is a product of small prime numbers.
!
!     wshigs an array which must be initialized by subroutine shigs .
!            once initialized, wshigs can be used repeatedly by shigs
!            as long as nlat and nlon remain unchanged.  wshigs must
!            not be altered between calls of shigs.
!
!     lshigs the dimension of the array wshigs as it appears in the
!            program that calls shigs. define
!
!               l1 = min(nlat, (nlon+2)/2) if nlon is even or
!               l1 = min(nlat, (nlon+1)/2) if nlon is odd
!
!            and
!
!               l2 = nlat/2        if nlat is even or
!               l2 = (nlat+1)/2    if nlat is odd
!
!            then lshigs must be at least
!
!            nlat*(3*(l1+l2)-2)+(l1-1)*(l2*(2*nlat-l1)-3*l1)/2+nlon+15
!
!     work   a real work space which need not be saved
!
!     lwork  the dimension of the array work as it appears in the
!            program that calls shigs. lwork must be at least
!            4*nlat*(nlat+2)+2 in the routine calling shigs
!
!     dwork   a real work array that does not have to be saved.
!
!     ldwork  the length of dwork in the calling routine.  ldwork must
!             be at least nlat*(nlat+4)
!
!     output parameter
!
!     wshags an array which must be initialized before calling shags or
!            once initialized, wshags can be used repeatedly by shags or
!            as long as nlat and nlon remain unchanged.  wshags must not
!            altered between calls of shasc.
!
!     ierror = 0  no errors
!            = 1  error in the specification of nlat
!            = 2  error in the specification of nlon
!            = 3  error in the specification of lshags
!            = 4  error in the specification of lwork
!            = 5  error in the specification of ldwork
!            = 6  failure in gaqd to compute gaussian points
!                 (due to failure in eigenvalue routine)
!
!
! ****************************************************************
!
subroutine shigs(nlat, nlon, wshigs, lshigs, work, lwork, dwork, &
    ldwork, ierror)
    !
    !     this subroutine must be called before calling shags or shsgs with
    !     fixed nlat, nlon. it precomputes the gaussian weights, points
    !     and all necessary legendre polys and stores them in wshigs.
    !     these quantities must be preserved when calling shsgs or shags
    !     repeatedly with fixed nlat, nlon.
    !
    dimension wshigs(lshigs), work(lwork)
    real dwork(ldwork)
    ierror = 1
    if (nlat<3) return
    ierror = 2
    if (nlon<4) return
    !     set triangular truncation limit for spherical harmonic basis
    l = min((nlon+2)/2, nlat)
    !     set equator or nearest point (if excluded) pointer
    late = (nlat+1)/2
    l1 = l
    l2 = late
    !     check permanent work space length
    ierror = 3
    lp=nlat*(3*(l1+l2)-2)+(l1-1)*(l2*(2*nlat-l1)-3*l1)/2+nlon+15
    if (lshigs<lp) return
    ierror = 4
    !     check temporary work space
    if (lwork<4*nlat*(nlat+2)+2) return
    !     check temp real space
    ierror = 5
    if (ldwork < nlat*(nlat+4)) return
    ierror = 0
    !     set preliminary quantites needed to compute and store legendre polys
    call shigsp(nlat, nlon, wshigs, lshigs, dwork, ldwork, ierror)
    if (ierror/=0) return
    !     set legendre poly pointer in wshigs
    ipmnf = nlat+2*nlat*late+3*(l*(l-1)/2+(nlat-l)*(l-1))+nlon+16
    call shigss1(nlat, l, late, wshigs, work, wshigs(ipmnf))
    return
end subroutine shigs

subroutine shigss1(nlat, l, late, w, pmn, pmnf)
    dimension w(1), pmn(nlat, late, 3), pmnf(late, 1)
    !     compute and store legendre polys for i=1, ..., late, m=0, ..., l-1
    !     and n=m, ..., l-1

    ! Initialize
    pmn = 0.0

    do mp1=1, l
        m = mp1-1
        mml1 = m*(2*nlat-m-1)/2
        !     compute pmn for n=m, ..., nlat-1 and i=1, ..., (l+1)/2
        mode = 0
        call legin(mode, l, nlat, m, w, pmn, km)
        !     store above in pmnf
        do np1=mp1, nlat
            mn = mml1+np1
            do i=1, late
                pmnf(i, mn) = pmn(np1, i, km)
            end do
        end do
    end do

end subroutine shigss1


subroutine shigsp(nlat, nlon, wshigs, lshigs, dwork, ldwork, ierror)
    dimension wshigs(lshigs)
    real dwork(ldwork)
    ierror = 1
    if (nlat<3) return
    ierror = 2
    if (nlon<4) return
    !     set triangular truncation limit for spherical harmonic basis
    l = min((nlon+2)/2, nlat)
    !     set equator or nearest point (if excluded) pointer
    late = (nlat+mod(nlat, 2))/2
    l1 = l
    l2 = late
    ierror = 3
    !     check permanent work space length
    if (lshigs < nlat*(2*l2+3*l1-2)+3*l1*(1-l1)/2+nlon+15)return
    ierror = 4
    !     if (lwork.lt.4*nlat*(nlat+2)+2) return
    if (ldwork < nlat*(nlat+4)) return
    ierror = 0
    !     set pointers
    i1 = 1
    i2 = i1+nlat
    i3 = i2+nlat*late
    i4 = i3+nlat*late
    i5 = i4+l*(l-1)/2 +(nlat-l)*(l-1)
    i6 = i5+l*(l-1)/2 +(nlat-l)*(l-1)
    i7 = i6+l*(l-1)/2 +(nlat-l)*(l-1)
    !     set indices in temp work for real gaussian wts and pts
    idth = 1
    !     idwts = idth+2*nlat
    !     iw = idwts+2*nlat
    idwts = idth+nlat
    iw = idwts+nlat
    call shigsp1(nlat, nlon, l, late, wshigs(i1), wshigs(i2), wshigs(i3), &
        wshigs(i4), wshigs(i5), wshigs(i6), wshigs(i7), dwork(idth), &
        dwork(idwts), dwork(iw), ierror)

    if (ierror/=0) then
        ierror = 5
    end if

end subroutine shigsp



subroutine shigsp1(nlat, nlon, l, late, wts, p0n, p1n, abel, bbel, cbel, &
    wfft, dtheta, dwts, work, ier)
    dimension wts(nlat), p0n(nlat, late), p1n(nlat, late), abel(1), bbel(1), &
        cbel(1), wfft(1), dtheta(nlat), dwts(nlat)
    real pb, dtheta, work(*)

    indx(m, n) = (n-1)*(n-2)/2+m-1
    imndx(m, n) = l*(l-1)/2+(n-l-1)*(l-1)+m-1

    call hrffti(nlon, wfft)
    !     compute real gaussian points and weights
    !     lw = 4*nlat*(nlat+2)
    lw = nlat*(nlat+2)

    call gaqd(nlat, dtheta, dwts, work, lw, ier)

    if (ier/=0) then
        return
    end if

    !     store gaussian weights single precision to save computation
    !     in inner loops in analysis
    wts = dwts

    ! initialize p0n, p1n using real dnlfk, dnlft
    p0n = 0.0
    p1n = 0.0

    !     compute m=n=0 legendre polynomials for all theta(i)
    np1 = 1
    n = 0
    m = 0

    call dnlfk(m, n, work)

    do i=1, late
        call dnlft(m, n, dtheta(i), work, pb)
        p0n(1, i) = pb
    end do

    !     compute p0n, p1n for all theta(i) when n.gt.0
    do np1=2, nlat

        n = np1-1
        m = 0

        call dnlfk(m, n, work)

        do i=1, late
            call dnlft(m, n, dtheta(i), work, pb)
            p0n(np1, i) = pb
        end do
        !     compute m=1 legendre polynomials for all n and theta(i)
        m = 1
        call dnlfk(m, n, work)
        do i=1, late
            call dnlft(m, n, dtheta(i), work, pb)
            p1n(np1, i) = pb
        end do
    end do
    !
    !     compute and store swarztrauber recursion coefficients
    !     for 2.le.m.le.n and 2.le.n.le.nlat in abel, bbel, cbel
    do n=2, nlat
        mlim = min(n, l)
        do m=2, mlim
            imn = indx(m, n)
            if (n >= l) imn = imndx(m, n)
            abel(imn)=sqrt(real((2*n+1)*(m+n-2)*(m+n-3))/ &
                real(((2*n-3)*(m+n-1)*(m+n))))
            bbel(imn)=sqrt(real((2*n+1)*(n-m-1)*(n-m))/ &
                real(((2*n-3)*(m+n-1)*(m+n))))
            cbel(imn)=sqrt(real((n-m+1)*(n-m+2))/ &
                real(((n+m-1)*(n+m))))
        end do
    end do

end subroutine shigsp1
