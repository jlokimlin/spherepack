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
subroutine lfin (init, theta, l, m, nm, pb, id, wlfin)
dimension       pb(1)        , wlfin(1)
!
!     total length of wlfin is 4*l*(nm+1)
!
lnx = l*(nm+1)
iw1 = lnx+1
iw2 = iw1+lnx
iw3 = iw2+lnx
call lfin1(init, theta, l, m, nm, id, pb, wlfin, wlfin(iw1), &
                wlfin(iw2), wlfin(iw3), wlfin(iw2))
return
end subroutine lfin
subroutine lfin1(init, theta, l, m, nm, id, p3, phz, ph1, p1, p2, cp)
dimension       p1(l, 1)    , p2(l, 1)    , p3(id, 1)   , phz(l, 1)   , &
                ph1(l, 1)   , cp(1)      , theta(1)
nmp1 = nm+1
if(init /= 0) go to 5
ssqrt2 = 1./sqrt(2.)
do 10 i=1, l
phz(i, 1) = ssqrt2
10 continue
do 15 np1=2, nmp1
nh = np1-1
call alfk(nh, 0, cp)
do 16 i=1, l
call lfpt(nh, 0, theta(i), cp, phz(i, np1))
16 continue
call alfk(nh, 1, cp)
do 17 i=1, l
call lfpt(nh, 1, theta(i), cp, ph1(i, np1))
17 continue
15 continue
return
5 mp1 = m+1
fm = real(m)
tm = fm+fm
if(m-1)25, 30, 35
25 do 45 np1=1, nmp1
do 45 i=1, l
p3(i, np1) = phz(i, np1)
p1(i, np1) = phz(i, np1)
45 continue
return
30 do 50 np1=2, nmp1
do 50 i=1, l
p3(i, np1) = ph1(i, np1)
p2(i, np1) = ph1(i, np1)
50 continue
return
35 temp = tm*(tm-1.)
cc = sqrt((tm+1.)*(tm-2.)/temp)
ee = sqrt(2./temp)
do 85 i=1, l
p3(i, m+1) = cc*p1(i, m-1)-ee*p1(i, m+1)
85 continue
if(m == nm) return
temp = tm*(tm+1.)
cc = sqrt((tm+3.)*(tm-2.)/temp)
ee = sqrt(6./temp)
do 70 i=1, l
p3(i, m+2) = cc*p1(i, m)-ee*p1(i, m+2)
70 continue
mp3 = m+3
if(nmp1 < mp3) go to 80
do 75 np1=mp3, nmp1
n = np1-1
fn = real(n)
tn = fn+fn
cn = (tn+1.)/(tn-3.)
fnpm = fn+fm
fnmm = fn-fm
temp = fnpm*(fnpm-1.)
cc = sqrt(cn*(fnpm-3.)*(fnpm-2.)/temp)
dd = sqrt(cn*fnmm*(fnmm-1.)/temp)
ee = sqrt((fnmm+1.)*(fnmm+2.)/temp)
do 75 i=1, l
p3(i, np1) = cc*p1(i, np1-2)+dd*p3(i, np1-2)-ee*p1(i, np1)
75 continue
80 do 90 np1=m, nmp1
do 90 i=1, l
p1(i, np1) = p2(i, np1)
p2(i, np1) = p3(i, np1)
90 continue
return
end subroutine lfin1

