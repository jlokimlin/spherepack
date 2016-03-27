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
! ... file shaes.f
!
!     this file contains code and documentation for subroutines
!     shaes and shaesi
!
! ... files which must be loaded with shaes.f
!
!     sphcom.f, hrfft.f
!
!     subroutine shaes(nlat, nlon, isym, nt, g, idg, jdg, a, b, mdab, ndab, 
!    +                 wshaes, lshaes, work, lwork, ierror)
!
!     subroutine shaes performs the spherical harmonic analysis
!     on the array g and stores the result in the arrays a and b.
!     the analysis is performed on an equally spaced grid.  the
!     associated legendre functions are stored rather than recomputed
!     as they are in subroutine shaec.  the analysis is described
!     below at output parameters a, b.
!
!     sphcom.f, hrfft.f
!
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
!            than or equal to 4. the efficiency of the computation is
!            improved when nlon is a product of small prime numbers.
!
!     isym   = 0  no symmetries exist about the equator. the analysis
!                 is performed on the entire sphere.  i.e. on the
!                 array g(i, j) for i=1, ..., nlat and j=1, ..., nlon.
!                 (see description of g below)
!
!            = 1  g is antisymmetric about the equator. the analysis
!                 is performed on the northern hemisphere only.  i.e.
!                 if nlat is odd the analysis is performed on the
!                 array g(i, j) for i=1, ..., (nlat+1)/2 and j=1, ..., nlon.
!                 if nlat is even the analysis is performed on the
!                 array g(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
!
!
!            = 2  g is symmetric about the equator. the analysis is
!                 performed on the northern hemisphere only.  i.e.
!                 if nlat is odd the analysis is performed on the
!                 array g(i, j) for i=1, ..., (nlat+1)/2 and j=1, ..., nlon.
!                 if nlat is even the analysis is performed on the
!                 array g(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
!
!     nt     the number of analyses.  in the program that calls shaes, 
!            the arrays g, a and b can be three dimensional in which
!            case multiple analyses will be performed.  the third
!            index is the analysis index which assumes the values
!            k=1, ..., nt.  for a single analysis set nt=1. the
!            discription of the remaining parameters is simplified
!            by assuming that nt=1 or that the arrays g, a and b
!            have only two dimensions.
!
!     g      a two or three dimensional array (see input parameter
!            nt) that contains the discrete function to be analyzed.
!            g(i, j) contains the value of the function at the colatitude
!            point theta(i) = (i-1)*pi/(nlat-1) and longitude point
!            phi(j) = (j-1)*2*pi/nlon. the index ranges are defined
!            above at the input parameter isym.
!
!
!     idg    the first dimension of the array g as it appears in the
!            program that calls shaes.  if isym equals zero then idg
!            must be at least nlat.  if isym is nonzero then idg
!            must be at least nlat/2 if nlat is even or at least
!            (nlat+1)/2 if nlat is odd.
!
!     jdg    the second dimension of the array g as it appears in the
!            program that calls shaes.  jdg must be at least nlon.
!
!     mdab   the first dimension of the arrays a and b as it appears
!            in the program that calls shaes. mdab must be at least
!            min(nlat, (nlon+2)/2) if nlon is even or at least
!            min(nlat, (nlon+1)/2) if nlon is odd.
!
!     ndab   the second dimension of the arrays a and b as it appears
!            in the program that calls shaes. ndab must be at least nlat
!
!     wshaes an array which must be initialized by subroutine shaesi.
!            once initialized, wshaes can be used repeatedly by shaes
!            as long as nlon and nlat remain unchanged.  wshaes must
!            not be altered between calls of shaes.
!
!     lshaes the dimension of the array wshaes as it appears in the
!            program that calls shaes. define
!
!               l1 = min(nlat, (nlon+2)/2) if nlon is even or
!               l1 = min(nlat, (nlon+1)/2) if nlon is odd
!
!            and
!
!               l2 = nlat/2        if nlat is even or
!               l2 = (nlat+1)/2    if nlat is odd
!
!            then lshaes must be at least
!
!               (l1*l2*(nlat+nlat-l1+1))/2+nlon+15
!
!     work   a work array that does not have to be saved.
!
!     lwork  the dimension of the array work as it appears in the
!            program that calls shaes.  define
!
!               l2 = nlat/2        if nlat is even or
!               l2 = (nlat+1)/2    if nlat is odd
!
!            if isym is zero then lwork must be at least
!            (nt+1)*nlat*nlon. if isym is not zero then
!            lwork must be at least (nt+1)*l2*nlon.
!
!
!     **************************************************************
!
!     output parameters
!
!     a, b    both a, b are two or three dimensional arrays (see input
!            parameter nt) that contain the spherical harmonic
!            coefficients in the representation of g(i, j) given in the
!            discription of subroutine shses. for isym=0, a(m, n) and
!            b(m, n) are given by the equations listed below. symmetric
!            versions are used when isym is greater than zero.
!
!
!
!     definitions
!
!     1. the normalized associated legendre functions
!
!     pbar(m, n, theta) = sqrt((2*n+1)*factorial(n-m)/(2*factorial(n+m)))
!                       *sin(theta)**m/(2**n*factorial(n)) times the
!                       (n+m)th derivative of (x**2-1)**n with respect
!                       to x=cos(theta)
!
!     2. the normalized z functions for m even
!
!     zbar(m, n, theta) = 2/(nlat-1) times the sum from k=0 to k=nlat-1 of
!                       the integral from tau = 0 to tau = pi of
!                       cos(k*theta)*cos(k*tau)*pbar(m, n, tau)*sin(tau)
!                       (first and last terms in this sum are divided
!                       by 2)
!
!     3. the normalized z functions for m odd
!
!     zbar(m, n, theta) = 2/(nlat-1) times the sum from k=0 to k=nlat-1 of
!                       of the integral from tau = 0 to tau = pi of
!                       sin(k*theta)*sin(k*tau)*pbar(m, n, tau)*sin(tau)
!
!     4. the fourier transform of g(i, j).
!
!     c(m, i)          = 2/nlon times the sum from j=1 to j=nlon
!                       of g(i, j)*cos((m-1)*(j-1)*2*pi/nlon)
!                       (the first and last terms in this sum
!                       are divided by 2)
!
!     s(m, i)          = 2/nlon times the sum from j=2 to j=nlon
!                       of g(i, j)*sin((m-1)*(j-1)*2*pi/nlon)
!
!     5. the maximum (plus one) longitudinal wave number
!
!            mmax = min(nlat, (nlon+2)/2) if nlon is even or
!            mmax = min(nlat, (nlon+1)/2) if nlon is odd.
!
!     then for m=0, ..., mmax-1 and  n=m, ..., nlat-1  the arrays a, b are
!     given by
!
!     a(m+1, n+1)      = the sum from i=1 to i=nlat of
!                       c(m+1, i)*zbar(m, n, theta(i))
!                       (first and last terms in this sum are
!                       divided by 2)
!
!     b(m+1, n+1)      = the sum from i=1 to i=nlat of
!                       s(m+1, i)*zbar(m, n, theta(i))
!
!
!     ierror = 0  no errors
!            = 1  error in the specification of nlat
!            = 2  error in the specification of nlon
!            = 3  error in the specification of isym
!            = 4  error in the specification of nt
!            = 5  error in the specification of idg
!            = 6  error in the specification of jdg
!            = 7  error in the specification of mdab
!            = 8  error in the specification of ndab
!            = 9  error in the specification of lshaes
!            = 10 error in the specification of lwork
!
!
! ****************************************************************
!     subroutine shaesi(nlat, nlon, wshaes, lshaes, work, lwork, dwork, 
!    +                  ldwork, ierror)
!
!     subroutine shaesi initializes the array wshaes which can then
!     be used repeatedly by subroutine shaes
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
!            than or equal to 4. the efficiency of the computation is
!            improved when nlon is a product of small prime numbers.
!
!     lshaes the dimension of the array wshaes as it appears in the
!            program that calls shaesi. define
!
!               l1 = min(nlat, (nlon+2)/2) if nlon is even or
!               l1 = min(nlat, (nlon+1)/2) if nlon is odd
!
!            and
!
!               l2 = nlat/2        if nlat is even or
!               l2 = (nlat+1)/2    if nlat is odd
!
!            then lshaes must be at least
!
!               (l1*l2*(nlat+nlat-l1+1))/2+nlon+15
!
!     work   a real   work array that does not have to be saved.
!
!     lwork  the dimension of the array work as it appears in the
!            program that calls shaesi.  define
!
!               l1 = min(nlat, (nlon+2)/2) if nlon is even or
!               l1 = min(nlat, (nlon+1)/2) if nlon is odd
!
!            and
!
!               l2 = nlat/2        if nlat is even or
!               l2 = (nlat+1)/2    if nlat is odd
!
!            then lwork must be at least
!
!               5*nlat*l2+3*((l1-2)*(nlat+nlat-l1-1))/2
!
!
!     dwork  a real work array that does not have to be saved.
!
!     ldwork the dimension of the array dwork as it appears in the
!            program that calls shaesi.  ldwork must be at least nlat+1
!
!
!     output parameters
!
!     wshaes an array which is initialized for use by subroutine shaes.
!            once initialized, wshaes can be used repeatedly by shaes
!            as long as nlon and nlat remain unchanged.  wshaes must
!            not be altered between calls of shaes.
!
!     ierror = 0  no errors
!            = 1  error in the specification of nlat
!            = 2  error in the specification of nlon
!            = 3  error in the specification of lshaes
!            = 4  error in the specification of lwork
!            = 5  error in the specification of ldwork
!
!
! ****************************************************************
subroutine shaes(nlat, nlon, isym, nt, g, idg, jdg, a, b, mdab, ndab, &
                    wshaes, lshaes, work, lwork, ierror)
dimension g(idg, jdg, 1), a(mdab, ndab, 1), b(mdab, ndab, 1), wshaes(1), &
          work(1)
ierror = 1
if(nlat<3) return
ierror = 2
if(nlon<4) return
ierror = 3
if(isym<0 .or. isym>2) return
ierror = 4
if(nt < 0) return
ierror = 5
if((isym==0 .and. idg<nlat) .or. &
   (isym/=0 .and. idg<(nlat+1)/2)) return
ierror = 6
if(jdg < nlon) return
ierror = 7
mmax = min(nlat, nlon/2+1)
if(mdab < mmax) return
ierror = 8
if(ndab < nlat) return
ierror = 9
imid = (nlat+1)/2
idz = (mmax*(nlat+nlat-mmax+1))/2
lzimn = idz*imid
if(lshaes < lzimn+nlon+15) return
ierror = 10
ls = nlat
if(isym > 0) ls = imid
nln = nt*ls*nlon
if(lwork < nln+ls*nlon) return
ierror = 0
ist = 0
if(isym == 0) ist = imid
call shaes1(nlat, isym, nt, g, idg, jdg, a, b, mdab, ndab, wshaes, idz, &
         ls, nlon, work, work(ist+1), work(nln+1), wshaes(lzimn+1))
return
end subroutine shaes
subroutine shaes1(nlat, isym, nt, g, idgs, jdgs, a, b, mdab, ndab, z, idz, &
                  idg, jdg, ge, go, work, whrfft)
dimension g(idgs, jdgs, 1), a(mdab, ndab, 1), b(mdab, ndab, 1), z(idz, 1), &
          ge(idg, jdg, 1), go(idg, jdg, 1), work(1), whrfft(1)
ls = idg
nlon = jdg
mmax = min(nlat, nlon/2+1)
mdo = mmax
if(mdo+mdo-1 > nlon) mdo = mmax-1
nlp1 = nlat+1
tsn = 2./nlon
fsn = 4./nlon
imid = (nlat+1)/2
modl = mod(nlat, 2)
imm1 = imid
if(modl /= 0) imm1 = imid-1
if(isym /= 0) go to 15
do 5 k=1, nt
do 5 i=1, imm1
do 5 j=1, nlon
ge(i, j, k) = tsn*(g(i, j, k)+g(nlp1-i, j, k))
go(i, j, k) = tsn*(g(i, j, k)-g(nlp1-i, j, k))
5 continue
go to 30
15 do 20 k=1, nt
do 20 i=1, imm1
do 20 j=1, nlon
ge(i, j, k) = fsn*g(i, j, k)
20 continue
if(isym == 1) go to 27
30 if(modl == 0) go to 27
do 25 k=1, nt
do 25 j=1, nlon
ge(imid, j, k) = tsn*g(imid, j, k)
25 continue
27 do 35 k=1, nt
call hrfftf(ls, nlon, ge(1, 1, k), ls, whrfft, work)
if(mod(nlon, 2) /= 0) go to 35
do 36 i=1, ls
ge(i, nlon, k) = .5*ge(i, nlon, k)
36 continue
35 continue
do 40 k=1, nt
do 40 mp1=1, mmax
do 40 np1=mp1, nlat
a(mp1, np1, k) = 0.
b(mp1, np1, k) = 0.
40 continue
if(isym == 1) go to 145
do 110 k=1, nt
do 110 i=1, imid
do 110 np1=1, nlat, 2
a(1, np1, k) = a(1, np1, k)+z(np1, i)*ge(i, 1, k)
110 continue
ndo = nlat
if(mod(nlat, 2) == 0) ndo = nlat-1
do 120 mp1=2, mdo
m = mp1-1
mb = m*(nlat-1)-(m*(m-1))/2
do 120 k=1, nt
do 120 i=1, imid
do 120 np1=mp1, ndo, 2
a(mp1, np1, k) = a(mp1, np1, k)+z(np1+mb, i)*ge(i, 2*mp1-2, k)
b(mp1, np1, k) = b(mp1, np1, k)+z(np1+mb, i)*ge(i, 2*mp1-1, k)
120 continue
if(mdo == mmax .or. mmax > ndo) go to 135
mb = mdo*(nlat-1)-(mdo*(mdo-1))/2
do 130 k=1, nt
do 130 i=1, imid
do 130 np1=mmax, ndo, 2
a(mmax, np1, k) = a(mmax, np1, k)+z(np1+mb, i)*ge(i, 2*mmax-2, k)
130 continue
135 if(isym == 2) return
145 do 150 k=1, nt
do 150 i=1, imm1
do 150 np1=2, nlat, 2
a(1, np1, k) = a(1, np1, k)+z(np1, i)*go(i, 1, k)
150 continue
ndo = nlat
if(mod(nlat, 2) /= 0) ndo = nlat-1
do 160 mp1=2, mdo
m = mp1-1
mp2 = mp1+1
mb = m*(nlat-1)-(m*(m-1))/2
do 160 k=1, nt
do 160 i=1, imm1
do 160 np1=mp2, ndo, 2
a(mp1, np1, k) = a(mp1, np1, k)+z(np1+mb, i)*go(i, 2*mp1-2, k)
b(mp1, np1, k) = b(mp1, np1, k)+z(np1+mb, i)*go(i, 2*mp1-1, k)
160 continue
mp2 = mmax+1
if(mdo == mmax .or. mp2 > ndo) return
mb = mdo*(nlat-1)-(mdo*(mdo-1))/2
do 170 k=1, nt
do 170 i=1, imm1
do 170 np1=mp2, ndo, 2
a(mmax, np1, k) = a(mmax, np1, k)+z(np1+mb, i)*go(i, 2*mmax-2, k)
170 continue
return
end subroutine shaes1

subroutine shaesi(nlat, nlon, wshaes, lshaes, work, lwork, dwork, &
                  ldwork, ierror)
dimension wshaes(*), work(*)
real dwork(*)
!
!     length of wshaes is (l*(l+1)*imid)/2+nlon+15
!     length of work is 5*l*imid + 3*((l-3)*l+2)/2
!
ierror = 1
if(nlat<3) return
ierror = 2
if(nlon<4) return
ierror = 3
mmax = min(nlat, nlon/2+1)
imid = (nlat+1)/2
lzimn = (imid*mmax*(nlat+nlat-mmax+1))/2
if(lshaes < lzimn+nlon+15) return
ierror = 4
labc = 3*((mmax-2)*(nlat+nlat-mmax-1))/2
if(lwork < 5*nlat*imid + labc) return
ierror = 5
if (ldwork < nlat+1) return
ierror = 0
iw1 = 3*nlat*imid+1
idz = (mmax*(nlat+nlat-mmax+1))/2
call sea1(nlat, nlon, imid, wshaes, idz, work, work(iw1), dwork)
call hrffti(nlon, wshaes(lzimn+1))
return
end subroutine shaesi
