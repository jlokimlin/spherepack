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
! ... file shses.f
!
!     this file contains code and documentation for subroutines
!     shses and shsesi
!
! ... files which must be loaded with shses.f
!
!     sphcom.f, hrfft.f
!
!     subroutine shses(nlat, nlon, isym, nt, g, idg, jdg, a, b, mdab, ndab, 
!    +                 wshses, lshses, work, lwork, ierror)
!
!     subroutine shses performs the spherical harmonic synthesis
!     on the arrays a and b and stores the result in the array g.
!     the synthesis is performed on an equally spaced grid.  the
!     associated legendre functions are stored rather than recomputed
!     as they are in subroutine shsec.  the synthesis is described
!     below at output parameter g.
!
! *** required files from spherepack2
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
!     isym   = 0  no symmetries exist about the equator. the synthesis
!                 is performed on the entire sphere.  i.e. on the
!                 array g(i, j) for i=1, ..., nlat and j=1, ..., nlon.
!                 (see description of g below)
!
!            = 1  g is antisymmetric about the equator. the synthesis
!                 is performed on the northern hemisphere only.  i.e.
!                 if nlat is odd the synthesis is performed on the
!                 array g(i, j) for i=1, ..., (nlat+1)/2 and j=1, ..., nlon.
!                 if nlat is even the synthesis is performed on the
!                 array g(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
!
!
!            = 2  g is symmetric about the equator. the synthesis is
!                 performed on the northern hemisphere only.  i.e.
!                 if nlat is odd the synthesis is performed on the
!                 array g(i, j) for i=1, ..., (nlat+1)/2 and j=1, ..., nlon.
!                 if nlat is even the synthesis is performed on the
!                 array g(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
!
!     nt     the number of syntheses.  in the program that calls shses, 
!            the arrays g, a and b can be three dimensional in which
!            case multiple syntheses will be performed.  the third
!            index is the synthesis index which assumes the values
!            k=1, ..., nt.  for a single synthesis set nt=1. the
!            discription of the remaining parameters is simplified
!            by assuming that nt=1 or that the arrays g, a and b
!            have only two dimensions.
!
!     idg    the first dimension of the array g as it appears in the
!            program that calls shses.  if isym equals zero then idg
!            must be at least nlat.  if isym is nonzero then idg
!            must be at least nlat/2 if nlat is even or at least
!            (nlat+1)/2 if nlat is odd.
!
!     jdg    the second dimension of the array g as it appears in the
!            program that calls shses.  jdg must be at least nlon.
!
!     a, b    two or three dimensional arrays (see the input parameter
!            nt) that contain the coefficients in the spherical harmonic
!            expansion of g(i, j) given below at the definition of the
!            output parameter g.  a(m, n) and b(m, n) are defined for
!            indices m=1, ..., mmax and n=m, ..., nlat where mmax is the
!            maximum (plus one) longitudinal wave number given by
!            mmax = min(nlat, (nlon+2)/2) if nlon is even or
!            mmax = min(nlat, (nlon+1)/2) if nlon is odd.
!
!     mdab   the first dimension of the arrays a and b as it appears
!            in the program that calls shses. mdab must be at least
!            min(nlat, (nlon+2)/2) if nlon is even or at least
!            min(nlat, (nlon+1)/2) if nlon is odd.
!
!     ndab   the second dimension of the arrays a and b as it appears
!            in the program that calls shses. ndab must be at least nlat
!
!     wshses an array which must be initialized by subroutine shsesi.
!            once initialized, wshses can be used repeatedly by shses
!            as long as nlon and nlat remain unchanged.  wshses must
!            not be altered between calls of shses.
!
!     lshses the dimension of the array wshses as it appears in the
!            program that calls shses. define
!
!               l1 = min(nlat, (nlon+2)/2) if nlon is even or
!               l1 = min(nlat, (nlon+1)/2) if nlon is odd
!
!            and
!
!               l2 = nlat/2        if nlat is even or
!               l2 = (nlat+1)/2    if nlat is odd
!
!            then lshses must be at least
!
!               (l1*l2*(nlat+nlat-l1+1))/2+nlon+15
!
!     work   a work array that does not have to be saved.
!
!     lwork  the dimension of the array work as it appears in the
!            program that calls shses.  define
!
!               l2 = nlat/2        if nlat is even or
!               l2 = (nlat+1)/2    if nlat is odd
!
!            if isym is zero then lwork must be at least
!
!               (nt+1)*nlat*nlon
!
!            if isym is nonzero lwork must be at least
!
!               (nt+1)*l2*nlon.
!
!     **************************************************************
!
!     output parameters
!
!     g      a two or three dimensional array (see input parameter
!            nt) that contains the spherical harmonic synthesis of
!            the arrays a and b at the colatitude point theta(i) =
!            (i-1)*pi/(nlat-1) and longitude point phi(j) =
!            (j-1)*2*pi/nlon. the index ranges are defined above at
!            at the input parameter isym.  for isym=0, g(i, j) is
!            given by the the equations listed below.  symmetric
!            versions are used when isym is greater than zero.
!
!     the normalized associated legendre functions are given by
!
!     pbar(m, n, theta) = sqrt((2*n+1)*factorial(n-m)/(2*factorial(n+m)))
!                       *sin(theta)**m/(2**n*factorial(n)) times the
!                       (n+m)th derivative of (x**2-1)**n with respect
!                       to x=cos(theta)
!
!     define the maximum (plus one) longitudinal wave number
!     as   mmax = min(nlat, (nlon+2)/2) if nlon is even or
!          mmax = min(nlat, (nlon+1)/2) if nlon is odd.
!
!     then g(i, j) = the sum from n=0 to n=nlat-1 of
!
!                   .5*pbar(0, n, theta(i))*a(1, n+1)
!
!              plus the sum from m=1 to m=mmax-1 of
!
!                   the sum from n=m to n=nlat-1 of
!
!              pbar(m, n, theta(i))*(a(m+1, n+1)*cos(m*phi(j))
!                                    -b(m+1, n+1)*sin(m*phi(j)))
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
!            = 9  error in the specification of lshses
!            = 10 error in the specification of lwork
!
!
! ****************************************************************
!     subroutine shsesi(nlat, nlon, wshses, lshses, work, lwork, dwork, 
!    +                  ldwork, ierror)
!
!     subroutine shsesi initializes the array wshses which can then
!     be used repeatedly by subroutine shses.
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
!     lshses the dimension of the array wshses as it appears in the
!            program that calls shsesi. define
!
!               l1 = min(nlat, (nlon+2)/2) if nlon is even or
!               l1 = min(nlat, (nlon+1)/2) if nlon is odd
!
!            and
!
!               l2 = nlat/2        if nlat is even or
!               l2 = (nlat+1)/2    if nlat is odd
!
!            then lshses must be at least
!
!               (l1*l2*(nlat+nlat-l1+1))/2+nlon+15
!
!     work   a real   work array that does not have to be saved.
!
!     lwork  the dimension of the array work as it appears in
!            the program that calls shsesi.  define
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
!            program that calls shsesi.  ldwork must be at least nlat+1
!
!
!     output parameters
!
!     wshses an array which is initialized for use by subroutine shses.
!            once initialized, wshses can be used repeatedly by shses
!            as long as nlon and nlat remain unchanged.  wshses must
!            not be altered between calls of shses.
!
!     ierror = 0  no errors
!            = 1  error in the specification of nlat
!            = 2  error in the specification of nlon
!            = 3  error in the specification of lshses
!            = 4  error in the specification of lwork
!            = 5  error in the specification of ldwork
!
! ****************************************************************
subroutine shses(nlat, nlon, isym, nt, g, idg, jdg, a, b, mdab, ndab, &
                    wshses, lshses, work, lwork, ierror)
dimension g(idg, jdg, 1), a(mdab, ndab, 1), b(mdab, ndab, 1), wshses(1), &
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
lpimn = (imid*mmax*(nlat+nlat-mmax+1))/2
if(lshses < lpimn+nlon+15) return
ierror = 10
ls = nlat
if(isym > 0) ls = imid
nln = nt*ls*nlon
if(lwork< nln+ls*nlon) return
ierror = 0
ist = 0
if(isym == 0) ist = imid
call shses1(nlat, isym, nt, g, idg, jdg, a, b, mdab, ndab, wshses, imid, &
       ls, nlon, work, work(ist+1), work(nln+1), wshses(lpimn+1))
return
end subroutine shses
subroutine shses1(nlat, isym, nt, g, idgs, jdgs, a, b, mdab, ndab, p, imid, &
                  idg, jdg, ge, go, work, whrfft)
dimension g(idgs, jdgs, 1), a(mdab, ndab, 1), b(mdab, ndab, 1), p(imid, 1), &
          ge(idg, jdg, 1), go(idg, jdg, 1), work(1), whrfft(1)
ls = idg
nlon = jdg
mmax = min(nlat, nlon/2+1)
mdo = mmax
if(mdo+mdo-1 > nlon) mdo = mmax-1
nlp1 = nlat+1
modl = mod(nlat, 2)
imm1 = imid
if(modl /= 0) imm1 = imid-1
do 80 k=1, nt
do 80 j=1, nlon
do 80 i=1, ls
ge(i, j, k) = 0.
8000 continue
800 continue
80 continue
if(isym == 1) go to 125
do 100 k=1, nt
do 100 np1=1, nlat, 2
do 100 i=1, imid
ge(i, 1, k)=ge(i, 1, k)+a(1, np1, k)*p(i, np1)
100 continue
ndo = nlat
if(mod(nlat, 2) == 0) ndo = nlat-1
do 110 mp1=2, mdo
m = mp1-1
mb = m*(nlat-1)-(m*(m-1))/2
do 110 np1=mp1, ndo, 2
mn = mb+np1
do 110 k=1, nt
do 110 i=1, imid
ge(i, 2*mp1-2, k) = ge(i, 2*mp1-2, k)+a(mp1, np1, k)*p(i, mn)
ge(i, 2*mp1-1, k) = ge(i, 2*mp1-1, k)+b(mp1, np1, k)*p(i, mn)
110 continue
if(mdo == mmax .or. mmax > ndo) go to 122
mb = mdo*(nlat-1)-(mdo*(mdo-1))/2
do 120 np1=mmax, ndo, 2
mn = mb+np1
do 120 k=1, nt
do 120 i=1, imid
ge(i, 2*mmax-2, k) = ge(i, 2*mmax-2, k)+a(mmax, np1, k)*p(i, mn)
120 continue
122 if(isym == 2) go to 155
125 do 140 k=1, nt
do 140 np1=2, nlat, 2
do 140 i=1, imm1
go(i, 1, k)=go(i, 1, k)+a(1, np1, k)*p(i, np1)
140 continue
ndo = nlat
if(mod(nlat, 2) /= 0) ndo = nlat-1
do 150 mp1=2, mdo
mp2 = mp1+1
m = mp1-1
mb = m*(nlat-1)-(m*(m-1))/2
do 150 np1=mp2, ndo, 2
mn = mb+np1
do 150 k=1, nt
do 150 i=1, imm1
go(i, 2*mp1-2, k) = go(i, 2*mp1-2, k)+a(mp1, np1, k)*p(i, mn)
go(i, 2*mp1-1, k) = go(i, 2*mp1-1, k)+b(mp1, np1, k)*p(i, mn)
150 continue
mp2 = mmax+1
if(mdo == mmax .or. mp2 > ndo) go to 155
mb = mdo*(nlat-1)-(mdo*(mdo-1))/2
do 152 np1=mp2, ndo, 2
mn = mb+np1
do 152 k=1, nt
do 152 i=1, imm1
go(i, 2*mmax-2, k) = go(i, 2*mmax-2, k)+a(mmax, np1, k)*p(i, mn)
152 continue
155 do 160 k=1, nt
if(mod(nlon, 2) /= 0) go to 157
do 156 i=1, ls
ge(i, nlon, k) = 2.*ge(i, nlon, k)
156 continue
157 call hrfftb(ls, nlon, ge(1, 1, k), ls, whrfft, work)
160 continue
if(isym /= 0) go to 180
do 170 k=1, nt
do 170 j=1, nlon
do 175 i=1, imm1
g(i, j, k) = .5*(ge(i, j, k)+go(i, j, k))
g(nlp1-i, j, k) = .5*(ge(i, j, k)-go(i, j, k))
175 continue
if(modl == 0) go to 170
g(imid, j, k) = .5*ge(imid, j, k)
170 continue
return
180 do 185 k=1, nt
do 185 i=1, imid
do 185 j=1, nlon
g(i, j, k) = .5*ge(i, j, k)
185 continue
return
end subroutine shses1

subroutine shsesi(nlat, nlon, wshses, lshses, work, lwork, dwork, &
                  ldwork, ierror)
dimension wshses(*), work(*)
real dwork(*)
ierror = 1
if(nlat<3) return
ierror = 2
if(nlon<4) return
ierror = 3
mmax = min(nlat, nlon/2+1)
imid = (nlat+1)/2
lpimn = (imid*mmax*(nlat+nlat-mmax+1))/2
if(lshses < lpimn+nlon+15) return
ierror = 4
labc = 3*((mmax-2)*(nlat+nlat-mmax-1))/2
if(lwork < 5*nlat*imid + labc) return
ierror = 5
if (ldwork < nlat+1) return
ierror = 0
iw1 = 3*nlat*imid+1
call SES1(nlat, nlon, imid, wshses, work, WORK(iw1), dwork)
call hrffti(nlon, wshses(lpimn+1))
return
end subroutine shsesi
