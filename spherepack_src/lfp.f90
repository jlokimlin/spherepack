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
! subroutine lfp (init, n, m, l, cp, pb, w)
!
! dimension of           cp((n/2)+1), pb(l), w(5*l+41)
! arguments
!
! purpose                routine lfp uses coefficients computed by
!                        routine alfk to calculate the single precision
!                        normalized associated legendre function pbar(n, 
!                        m, theta) at colatitudes theta=(i-1)*pi/(l-1), 
!                        i=1, ..., l. subroutine lfp evaluates pbar
!                        using one of the following trigonometric
!                        expansions
!
!                        1) for n even and m even, pbar(m, n, theta) =
!                           .5*cp(1) plus the sum from k=1 to k=n/2
!                           of cp(k)*cos(2*k*th)
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
! usage                  call lfp(init, n, m, l, cp, pb, w)
!
! arguments
!
! on input               init
!                          = 0 initialization only
!                          = 1 compute pbar(n, m, theta)
!
!                          lfp call with init = 0 initializes array w;
!                          no values of pbar(n, m, theta) are computed.
!                          init=0 should be used on the first call, or
!                          if l or w values differ from those in the
!                          previous call.
!
!                        n
!                          nonnegative integer, less than l, specifying
!                          the degree of pbar(n, m, theta)
!
!                        m
!                          is the order of pbar(n, m, theta). m can be
!                          any integer however pbar(n, m, theta) = 0
!                          if abs(m) is greater than n and
!                          pbar(n, m, theta) = (-1)**m*pbar(n, -m, theta)
!                          for negative m.
!
!                        l
!                          number of colatitudes theta=(i-1)*pi/(l-1)
!                          for i=1, ..., l where l is greater than 1.
!                          l must be an odd integer.
!
!                        cp
!                          single precision array of length (n/2)+1
!                          containing coefficients computed by routine
!                          alfk
!
!                        w
!                          a single precision work array with at
!                          least 5*l+41 locations
!
! on output              pb
!                          single precision array of length l containing
!                          pbar(n, m, theta), theta=(i-1)*pi/(l-1) for i=1
!                          , ..., l.
!
!                        w
!                          a single precision array containing values
!                          which must not be destroyed if the next call
!                          will have the same value of input parameter n
!
! special conditions     calls to routine lfp must be preceded by an
!                        appropriate call to routine alfk.
!
! precision              single
!
! algorithm              the trigonometric series formula used by
!                        routine lfp to calculate pbar(n, m, theta) for
!                        theta=(i-1)*pi/(l-1), i=1, ..., n, depends on
!                        m and n as follows:
!
!                           1) for n even and m even, the formula is
!                              .5*cp(1) plus the sum from k=1 to k=n/2
!                              of cp(k)*cos(2*k*theta)
!                           2) for n even and m odd. the formula is
!                              the sum from k=1 to k=n/2 of
!                              cp(k)*sin(2*k*theta)
!                           3) for n odd and m even, the formula is
!                              the sum from k=1 to k=(n+1)/2 of
!                              cp(k)*cos((2*k-1)*theta)
!                           4) for n odd and m odd, the formula is
!                              the sum from k=1 to k=(n+1)/2 of
!                              cp(k)*sin((2*k-1)*theta)
!
! accuracy               comparison between routines lfp and double
!                        precision dlfp on the cray1 indicates greater
!                        accuracy for smaller values of input parameter
!                        n.  agreement to 12 places was obtained for
!                        n=10 and to 11 places for n=100.
!
! timing                 time per call to routine lfp is dependent on
!                        the input parameters l and n.
!
subroutine lfp (init, n, m, l, cp, pb, w)
dimension       cp(1)       , pb(1)    , w(1)
!
do 10 i=1, l
pb(i) = 0.
10 continue
ma = abs(m)
if (ma > n) return
iw1 = l+l+12
iw2 = iw1+3*(l+1)/2+15
call lfp1(init, n, ma, l, cp, pb, w, w(iw1), w(iw2))
return
end subroutine lfp
subroutine lfp1(init, n, m, l, cp, p, wsave1, wsave2, wsave3)
dimension cp(*), p(*), wsave1(*), wsave2(*), wsave3(*)
save lc, lq, ls
if (init/=0) go to 41
lc=(l+1)/2
ls=lc-2
lq=lc-1
call sinti(ls, wsave1)
call costi(lc, wsave2)
call cosqi(lq, wsave3)
return
41 if (n< 0) then
    goto 10
else if (n == 0) then 
    goto 10
else 
    goto 40
end if
10 if (m< 0) then
    goto 20
else if (m == 0) then 
    goto 20
else 
    goto 40
end if
20 ssqrt2 = 1./sqrt(2.)
do  30 i=1, l
   p(i) = ssqrt2
30 continue
return
40 ls2 = (l+1)/2
lm1 = l-1
np1 = n+1
pi = acos( -1.0 )
dt = pi/lm1
nmod = mod(n, 2)
mmod = mod(m, 2)
if (nmod< 0) then
    goto 50
else if (nmod == 0) then 
    goto 50
else 
    goto 120
end if
50 if (mmod< 0) then
    goto 60
else if (mmod == 0) then 
    goto 60
else 
    goto 90
end if
60 kdp = n/2+1
do 70 i=1, kdp
p(i)=.5*cp(i)
70 continue
p(lc)=p(lc)+p(lc)
call cost(lc, p, wsave2)
do 80 i=1, lc
lmi=l-i
p(lmi+1)=p(i)
80 continue
go to 190
90 kdp=n/2
do 100 i=1, kdp
p(i+1)=.5*cp(i)
100 continue
p(ls+2)=0.
call sint(ls, p(2), wsave1)
do 110 i=1, ls
lmi=l-i
p(lmi)=-p(i+1)
110 continue
p(l)=0.
go to 190
120 kdp=(n+1)/2
    if (mmod< 0) then
        goto 140
    else if (mmod == 0) then 
        goto 140
    else 
        goto 160
    end if
140 do 130 i=1, kdp
p(i)=.25*cp(i)
130 continue
call cosqb(lq, p, wsave3)
do 150 i=1, lq
lmi=l-i
p(lmi+1)=-p(i)
150 continue
go to 190
160 do 180 i=1, kdp
p(i+1)=.25*cp(i)
180 continue
call sinqb(lq, p(2), wsave3)
do 170 i=1, lq
lmi=l-i
p(lmi)=p(i+1)
170 continue
p(l)=0.
190 return
end subroutine lfp1

