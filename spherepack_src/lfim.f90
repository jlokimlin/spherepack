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
subroutine lfim (init, theta, l, n, nm, pb, id, wlfim)
dimension       pb(1)        , wlfim(1)
!
!     total length of wlfim is 4*l*(nm+1)
!
lnx = l*(nm+1)
iw1 = lnx+1
iw2 = iw1+lnx
iw3 = iw2+lnx
call lfim1(init, theta, l, n, nm, id, pb, wlfim, wlfim(iw1), &
                wlfim(iw2), wlfim(iw3), wlfim(iw2))
return
end subroutine lfim
subroutine lfim1(init, theta, l, n, nm, id, p3, phz, ph1, p1, p2, cp)
dimension       p1(l, *)    , p2(l, *)    , p3(id, *)   , phz(l, *)   , &
                ph1(l, *)   , cp(*)      , theta(*)
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
5 if(n > 2) go to 60
if(n-1)25, 30, 35
25 do 45 i=1, l
p3(i, 1)=phz(i, 1)
45 continue
return
30 do 50 i=1, l
p3(i, 1) = phz(i, 2)
p3(i, 2) = ph1(i, 2)
50 continue
return
35 sq5s6 = sqrt(5./6.)
sq1s6 = sqrt(1./6.)
do 55 i=1, l
p3(i, 1) = phz(i, 3)
p3(i, 2) = ph1(i, 3)
p3(i, 3) = sq5s6*phz(i, 1)-sq1s6*p3(i, 1)
p1(i, 1) = phz(i, 2)
p1(i, 2) = ph1(i, 2)
p2(i, 1) = phz(i, 3)
p2(i, 2) = ph1(i, 3)
p2(i, 3) = p3(i, 3)
55 continue
return
60 nm1 = n-1
np1 = n+1
fn = real(n)
tn = fn+fn
cn = (tn+1.)/(tn-3.)
do 65 i=1, l
p3(i, 1) = phz(i, np1)
p3(i, 2) = ph1(i, np1)
65 continue
if(nm1 < 3) go to 71
do 70 mp1=3, nm1
m = mp1-1
fm = real(m)
fnpm = fn+fm
fnmm = fn-fm
temp = fnpm*(fnpm-1.)
cc = sqrt(cn*(fnpm-3.)*(fnpm-2.)/temp)
dd = sqrt(cn*fnmm*(fnmm-1.)/temp)
ee = sqrt((fnmm+1.)*(fnmm+2.)/temp)
do 70 i=1, l
p3(i, mp1) = cc*p1(i, mp1-2)+dd*p1(i, mp1)-ee*p3(i, mp1-2)
70 continue
71 fnpm = fn+fn-1.
temp = fnpm*(fnpm-1.)
cc = sqrt(cn*(fnpm-3.)*(fnpm-2.)/temp)
ee = sqrt(6./temp)
do 75 i=1, l
p3(i, n) = cc*p1(i, n-2)-ee*p3(i, n-2)
75 continue
fnpm = fn+fn
temp = fnpm*(fnpm-1.)
cc = sqrt(cn*(fnpm-3.)*(fnpm-2.)/temp)
ee = sqrt(2./temp)
do 80 i=1, l
p3(i, n+1) = cc*p1(i, n-1)-ee*p3(i, n-1)
80 continue
do 90 mp1=1, np1
do 90 i=1, l
p1(i, mp1) = p2(i, mp1)
p2(i, mp1) = p3(i, mp1)
90 continue
return
end subroutine lfim1
