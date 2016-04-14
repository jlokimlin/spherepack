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
!     file alf.f contains subroutines alfk, lfim, lfim1, lfin, lfin1, lfpt
!     for computing normalized associated legendre polynomials
!
! subroutine alfk (n, m, cp)
!
! dimension of           real cp(n/2 + 1)
! arguments
!
! purpose                routine alfk computes single precision fourier
!                        coefficients in the trigonometric series
!                        representation of the normalized associated
!                        legendre function pbar(n, m, theta) for use by
!                        routines lfp and lfpt in calculating single
!                        precision pbar(n, m, theta).
!
!                        first define the normalized associated
!                        legendre functions
!
!                        pbar(m, n, theta) = sqrt((2*n+1)*factorial(n-m)
!                        /(2*factorial(n+m)))*sin(theta)**m/(2**n*
!                        factorial(n)) times the (n+m)th derivative of
!                        (x**2-1)**n with respect to x=cos(theta)
!
!                        where theta is colatitude.
!
!                        then subroutine alfk computes the coefficients
!                        cp(k) in the following trigonometric
!                        expansion of pbar(m, n, theta).
!
!                        1) for n even and m even, pbar(m, n, theta) =
!                           .5*cp(1) plus the sum from k=1 to k=n/2
!                           of cp(k+1)*cos(2*k*th)
!
!                        2) for n even and m odd, pbar(m, n, theta) =
!                           the sum from k=1 to k=n/2 of
!                           cp(k)*sin(2*k*th)
!
!                        3) for n odd and m even, pbar(m, n, theta) =
!                           the sum from k=1 to k=(n+1)/2 of
!                           cp(k)*cos((2*k-1)*th)
!
!                        4) for n odd and m odd,  pbar(m, n, theta) =
!                           the sum from k=1 to k=(n+1)/2 of
!                           cp(k)*sin((2*k-1)*th)
!
!
! usage                  call alfk(n, m, cp)
!
! arguments
!
! on input               n
!                          nonnegative integer specifying the degree of
!                          pbar(n, m, theta)
!
!                        m
!                          is the order of pbar(n, m, theta). m can be
!                          any integer however cp is computed such that
!                          pbar(n, m, theta) = 0 if abs(m) is greater
!                          than n and pbar(n, m, theta) = (-1)**m*
!                          pbar(n, -m, theta) for negative m.
!
! on output              cp
!                          single precision array of length (n/2)+1
!                          which contains the fourier coefficients in
!                          the trigonometric series representation of
!                          pbar(n, m, theta)
!
!
! special conditions     none
!
! precision              single
!
! algorithm              the highest order coefficient is determined in
!                        closed form and the remainig coefficients are
!                        determined as the solution of a backward
!                        recurrence relation.
!
! accuracy               comparison between routines alfk and double
!                        precision dalfk on the cray1 indicates
!                        greater accuracy for smaller values
!                        of input parameter n.  agreement to 14
!                        places was obtained for n=10 and to 13
!                        places for n=100.
!
subroutine alfk (n, m, cp)
    implicit none
    !----------------------------------------------------------------------
    ! Dictionary: calling arguments
    !----------------------------------------------------------------------
    integer, intent (in)    :: n
    integer, intent (in)    :: m
    real,    intent (out)   :: cp(n/2+1)
    !----------------------------------------------------------------------
    ! Dictionary: local variables
    !----------------------------------------------------------------------
    integer         :: i, l, ma, nex, nmms2
    real, parameter :: sc10 = 1024.0
    real, parameter :: sc20 = sc10**2
    real, parameter :: sc40 = sc20**2
    real            :: a1, b1, c1, t1, t2
    real            :: fk, cp2, pm1
    real            :: fden, fnmh, fnum, fnnp1, fnmsq
    !----------------------------------------------------------------------

    cp(1) = 0.0
    ma = abs(m)

    if (ma > n) then
        return
    end if

    if (n-1 < 0) then
        go to 2
    else if (n-1 == 0) then
        go to 3
    else 
        go to 5
    end if

2   cp(1) = sqrt(2.0)
    return

3   if (ma /= 0) then
        go to 4
    end if

    cp(1) = sqrt(1.5)
    return

4   cp(1) = sqrt(0.75)

    if (m == -1) then
        cp(1) = -cp(1)
    end if

    return

5   if (mod(n+ma, 2) /= 0) then
        go to 10
    end if

    nmms2 = (n-ma)/2
    fnum = n+ma+1
    fnmh = n-ma+1
    pm1 = 1.0
    go to 15

10  nmms2 = (n-ma-1)/2
    fnum = n+ma+2
    fnmh = n-ma+2
    pm1 = -1.0

15  t1 = 1.0/sc20
    nex = 20
    fden = 2.0
    if (nmms2 < 1) then
        go to 20
    end if

    do  i=1, nmms2
        t1 = fnum*t1/fden
        if (t1 > sc20) then
            t1 = t1/sc40
            nex = nex+40
        end if
        fnum = fnum+2.0
        fden = fden+2.0
    end do

20  t1 = t1/2.0**(n-1-nex)

    if (mod(ma/2, 2) /= 0) then
        t1 = -t1
    end if

    t2 = 1.0
    if (ma == 0) then
        go to 26
    end if

    do  i=1, ma
        t2 = fnmh*t2/(fnmh+pm1)
        fnmh = fnmh+2.
    end do

26  cp2 = t1*sqrt((real(n)+0.5)*t2)
    fnnp1 = real(n*(n+1))
    fnmsq = fnnp1-2.0*(ma**2)
    l = (n+1)/2

    if (mod(n, 2) == 0 .and. mod(ma, 2) == 0) then
        l = l+1
    end if

    cp(l) = cp2

    if (m >= 0) then
        go to 29
    end if

    if (mod(ma, 2) /= 0) cp(l) = -cp(l)

29  if (l <= 1) then
        return
    end if

    fk = n
    a1 = (fk-2.0)*(fk-1.0)-fnnp1
    b1 = 2.0*(fk*fk-fnmsq)
    cp(l-1) = b1*cp(l)/a1

30  l = l-1

    if (l <= 1) then
        return
    end if

    fk = fk-2.0
    a1 = (fk-2.0)*(fk-1.0)-fnnp1
    b1 = -2.0*(fk*fk-fnmsq)
    c1 = (fk+1.0)*(fk+2.0)-fnnp1
    cp(l-1) = -(b1*cp(l)+c1*cp(l+1))/a1
    go to 30

end subroutine alfk

! subroutine lfim (init, theta, l, n, nm, pb, id, wlfim)
!
! dimension of           theta(l),  pb(id, nm+1),  wlfim(4*l*(nm+1))
! arguments
!
! purpose                given n and l, routine lfim calculates
!                        the normalized associated legendre functions
!                        pbar(n, m, theta) for m=0, ..., n and theta(i)
!                        for i=1, ..., l where
!
!                        pbar(m, n, theta) = sqrt((2*n+1)*factorial(n-m)
!                        /(2*factorial(n+m)))*sin(theta)**m/(2**n*
!                        factorial(n)) times the (n+m)th derivative of
!                        (x**2-1)**n with respect to x=cos(theta)
!
! usage                  call lfim (init, theta, l, n, nm, pb, id, wlfim)
!
! arguments
! on input               init
!                        = 0
!                            initialization only - using parameters
!                            l, nm and array theta, subroutine lfim
!                            initializes array wlfim for subsequent
!                            use in the computation of the associated
!                            legendre functions pb. initialization
!                            does not have to be repeated unless
!                            l, nm, or array theta are changed.
!                        = 1
!                            subroutine lfim uses the array wlfim that
!                            was computed with init = 0 to compute pb.
!
!                        theta
!                          an array that contains the colatitudes
!                          at which the associated legendre functions
!                          will be computed. the colatitudes must be
!                          specified in radians.
!
!                        l
!                          the length of the theta array. lfim is
!                          vectorized with vector length l.
!
!                        n
!                          nonnegative integer, less than nm, specifying
!                          degree of pbar(n, m, theta). subroutine lfim
!                          must be called starting with n=0. n must be
!                          incremented by one in subsequent calls and
!                          must not exceed nm.
!
!                        nm
!                          the maximum value of n and m
!
!                        id
!                          the first dimension of the two dimensional
!                          array pb as it appears in the program that
!                          calls lfim. (see output parameter pb)
!
!                        wlfim
!                          an array with length 4*l*(nm+1) which
!                          must be initialized by calling lfim
!                          with init=0 (see parameter init)  it
!                          must not be altered between calls to
!                          lfim.
!
!
! on output              pb
!                          a two dimensional array with first
!                          dimension id in the program that calls
!                          lfim. the second dimension of pb must
!                          be at least nm+1. starting with n=0
!                          lfim is called repeatedly with n being
!                          increased by one between calls. on each
!                          call, subroutine lfim computes
!                          = pbar(m, n, theta(i)) for m=0, ..., n and
!                          i=1, ...l.
!
!                        wlfim
!                          array containing values which must not
!                          be altered unless l, nm or the array theta
!                          are changed in which case lfim must be
!                          called with init=0 to reinitialize the
!                          wlfim array.
!
! special conditions     n must be increased by one between calls
!                        of lfim in which n is not zero.
!
! precision              single
!
!
! algorithm              routine lfim calculates pbar(n, m, theta) using
!                        a four term recurrence relation. (unpublished
!                        notes by paul n. swarztrauber)
!
subroutine lfim(init, theta, l, n, nm, pb, id, wlfim)
    implicit none
    !----------------------------------------------------------------------
    ! Dictionary: calling arguments
    !----------------------------------------------------------------------
    integer, intent (in)     :: init
    real,    intent (in)     :: theta(*)
    integer, intent (in)     :: l
    integer, intent (in)     :: n
    integer, intent (in)     :: nm
    real,    intent (in out) :: pb(1)
    integer, intent (in)     :: id
    real,    intent (in out) :: wlfim(1)
    !----------------------------------------------------------------------

    !
    !     total length of wlfim is 4*l*(nm+1)
    !
    associate(lnx => l*(nm+1) )
        associate( iw1=> lnx+1)
            associate( iw2 => iw1+lnx)
                associate( iw3 => iw2+lnx)
                    call lfim1(init, theta, l, n, nm, id, pb, wlfim, wlfim(iw1), &
                        wlfim(iw2), wlfim(iw3), wlfim(iw2))
                end associate
            end associate
        end associate
    end associate

end subroutine lfim


subroutine lfim1(init, theta, l, n, nm, id, p3, phz, ph1, p1, p2, cp)
    implicit none
    !----------------------------------------------------------------------
    ! Dictionary: calling arguments
    !----------------------------------------------------------------------
    integer, intent (in)     :: init
    real,    intent (in)     :: theta(*)
    integer, intent (in)     :: l
    integer, intent (in)     :: n
    integer, intent (in)     :: nm
    integer, intent (in)     :: id
    real,    intent (in out) :: p3(id,*)
    real,    intent (in out) :: phz(l,*)
    real,    intent (in out) :: ph1(l,*)
    real,    intent (in out) :: p1(l,*)
    real,    intent (in out) :: p2(l,*)
    real,    intent (in out) :: cp(*)
    !----------------------------------------------------------------------
    ! Dictionary: local variables
    !----------------------------------------------------------------------
    integer         :: i, m, nm1, nh, mp1, np1, mp3, nmp1
    real            :: cc, dd, ee, cn, fm, fn, fnmm, fnpm
    real            :: tm, tn, temp
    real, parameter :: SQRT2 = sqrt(2.0)
    real, parameter :: SQRT5 = sqrt(5.0)
    real, parameter :: SQRT6 = sqrt(6.0)
    real, parameter :: ONE_OVER_SQRT2 = 1.0/SQRT2
    real, parameter :: ONE_OVER_SQRT6 = 1.0/SQRT6
    real, parameter :: SQRT5_OVER_SQRT6 = SQRT5/SQRT6
    !----------------------------------------------------------------------

    nmp1 = nm+1

    if (init /= 0) then
        go to 5
    end if

    phz(:, 1) = ONE_OVER_SQRT2

    do np1=2, nmp1
        nh = np1-1
        call alfk(nh, 0, cp)
        do i=1, l
            call lfpt(nh, 0, theta(i), cp, phz(i, np1))
        end do
        call alfk(nh, 1, cp)
        do  i=1, l
            call lfpt(nh, 1, theta(i), cp, ph1(i, np1))
        end do
    end do


    return
5   if (n > 2) then
        go to 60
    end if

    if (n-1< 0) then
        go to 25
    else if (n-1 == 0) then
        go to 30
    else 
        go to 35
    end if

25  p3(:, 1) = phz(:, 1)
    return

30  p3(:, 1) = phz(:, 2)
    p3(:, 2) = ph1(:, 2)
    return

35  p3(:, 1) = phz(:, 3)
    p3(:, 2) = ph1(:, 3)
    p3(:, 3) = SQRT5_OVER_SQRT6 * phz(:, 1) - ONE_OVER_SQRT6 * p3(:, 1)
    p1(:, 1) = phz(:, 2)
    p1(:, 2) = ph1(:, 2)
    p2(:, 1) = phz(:, 3)
    p2(:, 2) = ph1(:, 3)
    p2(:, 3) = p3(:, 3)
    return
60  nm1 = n-1
    np1 = n+1
    fn = real(n)
    tn = fn+fn
    cn = (tn+1.0)/(tn-3.0)
    p3(:, 1) = phz(:, np1)
    p3(:, 2) = ph1(:, np1)

    if (nm1 < 3) then
        go to 71
    end if

    do mp1=3, nm1
        m = mp1-1
        fm = real(m)
        fnpm = fn+fm
        fnmm = fn-fm
        temp = fnpm*(fnpm-1.0)
        cc = sqrt(cn*(fnpm-3.0)*(fnpm-2.0)/temp)
        dd = sqrt(cn*fnmm*(fnmm-1.0)/temp)
        ee = sqrt((fnmm+1.0)*(fnmm+2.0)/temp)
        p3(:, mp1) = cc*p1(i, mp1-2)+dd*p1(:, mp1)-ee*p3(:, mp1-2)
    end do

71  fnpm = fn+fn-1.
    temp = fnpm*(fnpm-1.0)
    cc = sqrt(cn*(fnpm-3.0)*(fnpm-2.0)/temp)
    ee = sqrt(6.0/temp)

    p3(:, n) = cc*p1(:, n-2)-ee*p3(:, n-2)

    fnpm = fn+fn
    temp = fnpm*(fnpm-1.0)
    cc = sqrt(cn*(fnpm-3.0)*(fnpm-2.0)/temp)
    ee = sqrt(2.0/temp)
    p3(:, n+1) = cc*p1(:, n-1)-ee*p3(:, n-1)


    do mp1=1, np1
        p1(:, mp1) = p2(:, mp1)
        p2(:, mp1) = p3(:, mp1)
    end do

end subroutine lfim1
! subroutine lfin (init, theta, l, m, nm, pb, id, wlfin)
!
! dimension of           theta(l),  pb(id, nm+1),  wlfin(4*l*(nm+1))
! arguments
!
! purpose                given m and l, routine lfin calculates
!                        the normalized associated legendre functions
!                        pbar(n, m, theta) for n=m, ..., nm and theta(i)
!                        for i=1, ..., l where
!
!                        pbar(m, n, theta) = sqrt((2*n+1)*factorial(n-m)
!                        /(2*factorial(n+m)))*sin(theta)**m/(2**n*
!                        factorial(n)) times the (n+m)th derivative of
!                        (x**2-1)**n with respect to x=cos(theta)
!
! usage                  call lfin (init, theta, l, m, nm, pb, id, wlfin)
!
! arguments
! on input               init
!                        = 0
!                            initialization only - using parameters
!                            l, nm and the array theta, subroutine lfin
!                            initializes the array wlfin for subsequent
!                            use in the computation of the associated
!                            legendre functions pb. initialization does
!                            not have to be repeated unless l, nm or
!                            the array theta are changed.
!                        = 1
!                            subroutine lfin uses the array wlfin that
!                            was computed with init = 0 to compute pb
!
!                        theta
!                          an array that contains the colatitudes
!                          at which the associated legendre functions
!                          will be computed. the colatitudes must be
!                          specified in radians.
!
!                        l
!                          the length of the theta array. lfin is
!                          vectorized with vector length l.
!
!                        m
!                          nonnegative integer, less than nm, specifying
!                          degree of pbar(n, m, theta). subroutine lfin
!                          must be called starting with n=0. n must be
!                          incremented by one in subsequent calls and
!                          must not exceed nm.
!
!                        nm
!                          the maximum value of n and m
!
!                        id
!                          the first dimension of the two dimensional
!                          array pb as it appears in the program that
!                          calls lfin. (see output parameter pb)
!
!                        wlfin
!                          an array with length 4*l*(nm+1) which
!                          must be initialized by calling lfin
!                          with init=0 (see parameter init)  it
!                          must not be altered between calls to
!                          lfin.
!
!
! on output              pb
!                          a two dimensional array with first
!                          dimension id in the program that calls
!                          lfin. the second dimension of pb must
!                          be at least nm+1. starting with m=0
!                          lfin is called repeatedly with m being
!                          increased by one between calls. on each
!                          call, subroutine lfin computes pb(i, n+1)
!                          = pbar(m, n, theta(i)) for n=m, ..., nm and
!                          i=1, ...l.
!
!                        wlfin
!                          array containing values which must not
!                          be altered unless l, nm or the array theta
!                          are changed in which case lfin must be
!                          called with init=0 to reinitialize the
!                          wlfin array.
!
! special conditions     m must be increased by one between calls
!                        of lfin in which m is not zero.
!
! precision              single
!
! algorithm              routine lfin calculates pbar(n, m, theta) using
!                        a four term recurrence relation. (unpublished
!                        notes by paul n. swarztrauber)
!
subroutine lfin(init, theta, l, m, nm, pb, id, wlfin)
    implicit none
    !----------------------------------------------------------------------
    ! Dictionary: calling arguments
    !----------------------------------------------------------------------
    integer, intent (in)     :: init
    real,    intent (in)     :: theta(*)
    integer, intent (in)     :: l
    integer, intent (in)     :: m
    integer, intent (in)     :: nm
    real,    intent (in out) :: pb(1)
    integer, intent (in)     ::  id
    real,    intent (in out) :: wlfin(1)
    !----------------------------------------------------------------------

    !     total length of wlfin is 4*l*(nm+1)
    !
    associate( lnx => l*(nm+1) )
        associate( iw1 => lnx+1 )
            associate( iw2 => iw1+lnx )
                associate( iw3 => iw2+lnx )
                    call lfin1(init, theta, l, m, nm, id, pb, wlfin, wlfin(iw1), &
                        wlfin(iw2), wlfin(iw3), wlfin(iw2))
                end associate
            end associate
        end associate
    end associate

end subroutine lfin


subroutine lfin1(init, theta, l, m, nm, id, p3, phz, ph1, p1, p2, cp)
    implicit none
    !----------------------------------------------------------------------
    ! Dictionary: calling arguments
    !----------------------------------------------------------------------
    integer, intent (in) :: init
    real,    intent (in) :: theta(*)!theta(1)
    integer, intent (in) :: l
    integer, intent (in) :: m
    integer, intent (in) :: nm
    integer, intent (in) :: id
    real,    intent (in out) :: p3(id, 1)
    real,    intent (in out) :: phz(l,1)
    real,    intent (in out) :: ph1(l,1)
    real,    intent (in out) :: p1(l,1)
    real,    intent (in out) :: p2(l,1)
    real,    intent (in out) :: cp(1)
    !----------------------------------------------------------------------
    ! Dictionary: local variables
    !----------------------------------------------------------------------
    integer         :: i, n, nh, mp1, np1, mp3, nmp1
    real            :: cc, dd, ee, cn, fm, fn, fnmm, fnpm
    real            :: tm, tn, temp
    real, parameter :: SQRT2 = sqrt(2.0)
    real, parameter :: ONE_OVER_SQRT2 = 1.0/sqrt2
    !----------------------------------------------------------------------
    nmp1 = nm+1

    if (init /= 0) then
        go to 5
    end if

    phz(:, 1) = ONE_OVER_SQRT2

    do np1=2, nmp1
        nh = np1-1
        call alfk(nh, 0, cp)
        do i=1, l
            call lfpt(nh, 0, theta(i), cp, phz(i, np1))
        end do
        call alfk(nh, 1, cp)
        do i=1, l
            call lfpt(nh, 1, theta(i), cp, ph1(i, np1))
        end do
    end do
    return

5   mp1 = m+1
    fm = real(m)
    tm = fm+fm

    if (m-1 < 0) then
        go to 25
    else if (m-1 == 0) then
        go to 30
    else 
        go to 35
    end if

    25 do np1=1, nmp1
        p3(:, np1) = phz(:, np1)
        p1(:, np1) = phz(:, np1)
    end do
    return

    30 do np1=2, nmp1
        p3(:, np1) = ph1(:, np1)
        p2(:, np1) = ph1(:, np1)
    end do
    return

35  temp = tm*(tm-1.0)
    cc = sqrt((tm+1.0)*(tm-2.0)/temp)
    ee = sqrt(2.0/temp)
    p3(:, m+1) = cc*p1(:, m-1)-ee*p1(:, m+1)

    if (m == nm) then
        return
    end if

    temp = tm*(tm+1.0)
    cc = sqrt((tm+3.0)*(tm-2.0)/temp)
    ee = sqrt(6.0/temp)
    p3(:, m+2) = cc*p1(:, m)-ee*p1(:, m+2)
    mp3 = m+3

    if (nmp1 < mp3) then
        go to 80
    end if

    do np1=mp3, nmp1
        n = np1-1
        fn = real(n)
        tn = fn+fn
        cn = (tn+1.0)/(tn-3.0)
        fnpm = fn+fm
        fnmm = fn-fm
        temp = fnpm*(fnpm-1.0)
        cc = sqrt(cn*(fnpm-3.0)*(fnpm-2.0)/temp)
        dd = sqrt(cn*fnmm*(fnmm-1.0)/temp)
        ee = sqrt((fnmm+1.0)*(fnmm+2.0)/temp)
        p3(:, np1) = cc*p1(:, np1-2)+dd*p3(:, np1-2)-ee*p1(:, np1)
    end do

    80 do np1=m, nmp1
        p1(:, np1) = p2(:, np1)
        p2(:, np1) = p3(:, np1)
    end do

end subroutine lfin1


! subroutine lfpt (n, m, theta, cp, pb)
!
! dimension of
! arguments
!                        cp((n/2)+1)
!
! purpose                routine lfpt uses coefficients computed by
!                        routine alfk to compute the single precision
!                        normalized associated legendre function pbar(n, 
!                        m, theta) at colatitude theta.
!
! usage                  call lfpt(n, m, theta, cp, pb)
!
! arguments
!
! on input               n
!                          nonnegative integer specifying the degree of
!                          pbar(n, m, theta)
!                        m
!                          is the order of pbar(n, m, theta). m can be
!                          any integer however pbar(n, m, theta) = 0
!                          if abs(m) is greater than n and
!                          pbar(n, m, theta) = (-1)**m*pbar(n, -m, theta)
!                          for negative m.
!
!                        theta
!                          single precision colatitude in radians
!
!                        cp
!                          single precision array of length (n/2)+1
!                          containing coefficients computed by routine
!                          alfk
!
! on output              pb
!                          single precision variable containing
!                          pbar(n, m, theta)
!
! special conditions     calls to routine lfpt must be preceded by an
!                        appropriate call to routine alfk.
!
! precision              single
!
! algorithm              the trigonometric series formula used by
!                        routine lfpt to calculate pbar(n, m, th) at
!                        colatitude th depends on m and n as follows:
!
!                           1) for n even and m even, the formula is
!                              .5*cp(1) plus the sum from k=1 to k=n/2
!                              of cp(k)*cos(2*k*th)
!                           2) for n even and m odd. the formula is
!                              the sum from k=1 to k=n/2 of
!                              cp(k)*sin(2*k*th)
!                           3) for n odd and m even, the formula is
!                              the sum from k=1 to k=(n+1)/2 of
!                              cp(k)*cos((2*k-1)*th)
!                           4) for n odd and m odd, the formula is
!                              the sum from k=1 to k=(n+1)/2 of
!                              cp(k)*sin((2*k-1)*th)
!
! accuracy               comparison between routines lfpt and double
!                        precision dlfpt on the cray1 indicates greater
!                        accuracy for greater values on input parameter
!                        n.  agreement to 13 places was obtained for
!                        n=10 and to 12 places for n=100.
!
! timing                 time per call to routine lfpt is dependent on
!                        the input parameter n.
!
subroutine lfpt(n, m, theta, cp, pb)
    implicit none
    !----------------------------------------------------------------------
    ! Dictionary: calling arguments
    !----------------------------------------------------------------------
    integer, intent (in)     :: n
    integer, intent (in)     :: m
    real,    intent (in)     :: theta
    real,    intent (in out) :: cp(1)
    real,    intent (in out) :: pb
    !----------------------------------------------------------------------
    ! Dictionary: calling arguments
    !----------------------------------------------------------------------
    integer :: i,  k, ma, mmod, nmod, kp1, np1, kdo
    real    :: ct, st, cdt, cth, sdt, summation
    !----------------------------------------------------------------------

    pb = 0.0
    ma = abs(m)

    if (ma > n) then
        return
    end if

    if (n < 0) then
        go to 10
    else if (n == 0) then
        go to 10
    else 
        go to 30
    end if

10  if (ma < 0) then
        go to 20
    else if (ma == 0) then
        go to 20
    else 
        go to 30
    end if

20  pb= sqrt(0.5)
    go to 140

30  np1 = n+1
    nmod = mod(n, 2)
    mmod = mod(ma, 2)

    if (nmod < 0) then
        go to 40
    else if (nmod == 0) then
        go to 40
    else 
        go to 90
    end if

40  if (mmod < 0) then
        go to 50
    else if (mmod == 0) then
        go to 50
    else 
        go to 70
    end if

50  kdo = n/2+1
    cdt = cos(2.0*theta)
    sdt = sin(2.0*theta)
    ct = 1.0
    st = 0.0
    summation = 0.5*cp(1)

    do kp1=2, kdo
        cth = cdt*ct-sdt*st
        st = sdt*ct+cdt*st
        ct = cth
        summation = summation+cp(kp1)*ct
    end do

    pb= summation
    go to 140

70  kdo = n/2
    cdt = cos(2.0*theta)
    sdt = sin(2.0*theta)
    ct = 1.
    st = 0.
    summation = 0.0
    do  k=1, kdo
        cth = cdt*ct-sdt*st
        st = sdt*ct+cdt*st
        ct = cth
        summation = summation+cp(k)*st
    end do
    pb= summation
    go to 140

90  kdo = (n+1)/2

    if (mmod < 0) then
        go to 100
    else if (mmod == 0) then
        go to 100
    else 
        go to 120
    end if

100 cdt = cos(2.0*theta)
    sdt = sin(2.0*theta)
    ct = cos(theta)
    st = -sin(theta)
    summation = 0.0

    do k=1, kdo
        cth = cdt*ct-sdt*st
        st = sdt*ct+cdt*st
        ct = cth
        summation = summation+cp(k)*ct
    end do

    pb= summation
    go to 140

120 cdt = cos(2.0*theta)
    sdt = sin(2.0*theta)
    ct = cos(theta)
    st = -sin(theta)
    summation = 0.0

    do k=1, kdo
        cth = cdt*ct-sdt*st
        st = sdt*ct+cdt*st
        ct = cth
        summation = summation+cp(k)*st
    end do

    pb= summation

140 return

end subroutine lfpt
