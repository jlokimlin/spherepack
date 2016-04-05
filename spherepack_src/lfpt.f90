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
subroutine lfpt (n, m, theta, cp, pb)
dimension       cp(1)
!
pb = 0.
ma = abs(m)
if(ma > n) return
if(n< 0) then
    goto 10
else if(n == 0) then 
    goto 10
else 
    goto 30
end if
10 if(ma< 0) then
    goto 20
else if(ma == 0) then 
    goto 20
else 
    goto 30
end if
20 pb= sqrt(.5)
go to 140
30 np1 = n+1
nmod = mod(n, 2)
mmod = mod(ma, 2)
if(nmod< 0) then
    goto 40
else if(nmod == 0) then 
    goto 40
else 
    goto 90
end if
40 if(mmod< 0) then
    goto 50
else if(mmod == 0) then 
    goto 50
else 
    goto 70
end if
50 kdo = n/2+1
cdt = cos(theta+theta)
sdt = sin(theta+theta)
ct = 1.
st = 0.
sum = .5*cp(1)
do  60 kp1=2, kdo
   cth = cdt*ct-sdt*st
   st = sdt*ct+cdt*st
   ct = cth
   sum = sum+cp(kp1)*ct
60 continue
pb= sum
go to 140
70 kdo = n/2
cdt = cos(theta+theta)
sdt = sin(theta+theta)
ct = 1.
st = 0.
sum = 0.
do  80 k=1, kdo
   cth = cdt*ct-sdt*st
   st = sdt*ct+cdt*st
   ct = cth
   sum = sum+cp(k)*st
80 continue
pb= sum
go to 140
90 kdo = (n+1)/2
   if(mmod< 0) then
       goto 100
   else if(mmod == 0) then 
       goto 100
   else 
       goto 120
   end if
100 cdt = cos(theta+theta)
sdt = sin(theta+theta)
ct = cos(theta)
st = -sin(theta)
sum = 0.
do 110 k=1, kdo
   cth = cdt*ct-sdt*st
   st = sdt*ct+cdt*st
   ct = cth
   sum = sum+cp(k)*ct
110 continue
pb= sum
go to 140
120 cdt = cos(theta+theta)
sdt = sin(theta+theta)
ct = cos(theta)
st = -sin(theta)
sum = 0.
do 130 k=1, kdo
   cth = cdt*ct-sdt*st
   st = sdt*ct+cdt*st
   ct = cth
   sum = sum+cp(k)*st
130 continue
pb= sum
140 return
end subroutine lfpt

