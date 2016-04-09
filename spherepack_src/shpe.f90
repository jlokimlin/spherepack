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
!                            August 2003
!
! ... file shpe.f
!
!     this file contains code and documentation for subroutines
!     shpei and shpe.
!
! ... files which must be loaded with shpe.f
!
!     sphcom.f, hrfft.f
!
!     subroutine shpei initializes arrays wshp and iwshp for
!     subsequent repeated use by subroutine shpe, which
!     performs the harmonic projection equivalent to a
!     harmonic analysis followed by harmonic synthesis
!     but faster and with less memory. (see description of
!     subroutine shpe below)
!
!     subroutine shpei(nlat, nlon, isym, mtrunc, wshp, lwshp, iwshp, 
!    1 liwshp, work, lwork, ierror)
!
!     shpei initializes arrays wshp and iwshp for repeated use
!     by subroutine shpe ....
!
!     input parameters
!
!     nlat   the number of colatitudes on the full sphere including the
!            poles. for example, nlat = 37 for a five degree grid.
!            nlat determines the grid increment in colatitude as
!            pi/(nlat-1).  if nlat is odd the equator is located at
!            grid point i=(nlat+1)/2. if nlat is even the equator is
!            located half way between points i=nlat/2 and i=nlat/2+1.
!            nlat must be at least 3.
!
!     nlon   the number of distinct londitude points.  nlon determines
!            the grid increment in longitude as 2*pi/nlon. for example
!            nlon = 72 for a five degree grid. nlon must be greater
!            than or equal to 4. the efficiency of the computation is
!            improved when nlon is a product of small prime numbers.
!            nlon must beat least 4.
!
!     isym   currently not used.    
!
!     mtrunc the highest longitudinal wave number retained in the
!            projection. It must be less than or equal to
!            the minimum of nlat-1 and nlon/2. The first wave
!            number is zero. For example, if wave numbers 0 and
!            1 are desired then mtrunc = 1.
!
!     lwshp  the dimension of the array wshp as it appears in the
!            program that calls shpei. It must be at least
!            2*(nlat+1)**2+nlon+log2(nlon)
!
!     liwshp the dimension of the array iwshp as it appears in the
!            program that calls shpei. It must be at least
!            4*(nlat+1).
!
!     lwork  the dimension of the array work as it appears in the
!            program that calls shpei. It must be at least
!            1.25*(nlat+1)**2+7*nlat+8.
!
!     **************************************************************
!
!     output parameters
!
!     wshp   a single precision array that must be saved for
!            repeated use by subroutine shpe.        
!
!     iwshp  an integer array that must be saved for repeated
!            use by subroutine shpe.        
!
!     work   a real work array that does 
!            not have to be saved.
!
!     ierror = 0  no errors
!            = 1  error in the specification of nlat
!            = 2  error in the specification of nlon
!            = 3  error in the specification of isym
!            = 4  error in the specification of mtrunc
!            = 5  error in the specification of lwshp
!            = 6  error in the specification of liwshp
!            = 7  error in the specification of lwork
!
subroutine shpei(nlat, nlon, isym, mtrunc, wshp, lwshp, iwshp, &
 liwshp, work, lwork, ierror)
real work(*)
dimension wshp(*), iwshp(*)
!
ierror = 1
if (nlat < 3) return
ierror = 2
if (nlon < 4) return
!      ierror = 3
!      if (isym.lt.0 .or. isym.gt.2) return
ierror = 4
mmax = min(nlat-1, nlon/2)
if (mtrunc<0 .or. mtrunc>mmax) return
ierror = 5
lw1 = 2*(nlat+1)**2
log2n = log(real(nlon))/log(2.0)
if (lwshp<lw1+nlon+log2n) return
ierror = 6
if (liwshp<4*(nlat+1)) return
ierror = 7
mlwk = 1.25*(nlat+1)**2+7*nlat+8
if (lwork <mlwk) return
ierror = 0
!
call hrffti(nlon, wshp(lw1+1))
!
nte = (nlat+1)/2
nloc1 = 2*nte*nte
nloc2 = nlat+1
iw1 = 1
iw2 = iw1+nloc1
iw3 = iw2+nloc1
iw4 = iw3+nloc1
jw1 = 1
jw2 = jw1+nloc2
jw3 = jw2+nloc2
jw4 = jw3+nloc2
kw1 = 1
kw2 = kw1+nte
kw3 = kw2+nte
kw4 = kw3+nte
kw5 = kw4+nte+1
kw6 = kw5+nte
kw7 = kw6+nte
kw8 = kw7+nte
kw9 = kw8+nte
kw10 = kw9+nloc2+nloc2
kw11 = kw10+nloc2

kw12 = kw11+nloc1
kw13 = kw12+nloc1
!
call shpei1(nlat, nlon, isym, mtrunc, nte, ierror, wshp(iw1), wshp(iw2), &
  wshp(iw3), wshp(iw4), iwshp(jw1), iwshp(jw2), iwshp(jw3), &
  iwshp(jw4), work(kw1), work(kw2), work(kw3), work(kw4), work(kw5), &
  work(kw6), work(kw7), work(kw8), work(kw9), work(kw10), work(kw11), &
  work(kw12), work(kw11), work(kw12), work(kw13))
return
end subroutine shpei
subroutine shpei1(nlat, nlon, isym, mtrunc, idp, ierror, &
  pe, po, ze, zo, ipse, jzse, ipso, jzso, &
  cp, work, wx, s, e, thet, xx, z, a, b, we, ped, wo, pod, u)
!
real sum, eps, pi, dthet, v(1,1), a1, b1, c1
parameter (eps=epsilon(1.0))
real cp(idp), work(idp), wx(idp), s(idp+1), &
  e(idp), thet(idp), xx(idp), z(idp), u(idp, idp), &
  we(idp, idp, 2), ped(idp, idp, 2), a(4*idp), b(2*idp), &
  wo(idp, idp, 2), pod(idp, idp, 2)
!
dimension pe(idp, idp, 2), po(idp, idp, 2), ze(idp, idp, 2), &
  zo(idp, idp, 2), &
  ipse(idp, 2), jzse(idp, 2), ipso(idp, 2), jzso(idp, 2), &
  nshe(2), nsho(2)
!
ns2 = nlat/2
modn = nlat-ns2-ns2
nte = (nlat+1)/2
nto = nlat-nte
tusl = 0.
toe = 0.
!
!     compute grid distribution
!
pi = 4.0*atan(1.0)
dthet = pi/(nlat-1)
do i=1, nte
thet(i) = (i-1)*dthet
end do
!
!     compute weight matrices for even functions
!
do 40 mp1=1, 2
m = mp1-1
mrank = nlat-m-m
nem = (mrank+1)/2
do j=1, nem
n = j+j+m-2
call dlfkp(m, n, cp)
do i=1, nte
call dlftp (m, n, thet(i), cp, ped(i, j, mp1))
end do
if (m>0) ped(1, j, mp1) = 0.0
end do
call dsvdc(ped(m+1, 1, mp1), idp, nem, nem, s, e, u, &
                        idp, v(1,1), idp, work, 10, info)
!
do j=1, nem
s(j) = 1.0/(s(j)*s(j)) 
end do
!
!     compute weight matrix as u  s sup -2 u transpose
!
do j=1, nte
do i=1, nte
we(i, j, mp1) = 0.0
end do
end do
do i=1, nem
do j=1, nem
sum = 0.
do k=1, nem
sum = sum+s(k)*u(i, k)*u(j, k)
end do
we(i+m, j+m, mp1) = sum
end do
end do
40 continue
we(1, 1, 2) = 1.0
!
!     compute n**2 basis (even functions)
!
do n=1, nlat+nlat-2
dfn = n
a(n) = sqrt(dfn*(dfn+1.0))
end do
do n=1, nlat-1
dfn = n
b(n) = sqrt((dfn+dfn+3.0)/(dfn+dfn-1.0))
end do
!
mxtr = min(nlat-1, nlon/2, mtrunc)
ip = 2
do 200 mp1=1, mxtr+1
m = mp1-1
ip = 3-ip
ms2 = mp1/2
nrank = ms2+ms2
mrank = nlat-nrank
nem = (mrank+1)/2
!
!     compute associated legendre functions
!
if (m<=1) then
do 205 j=1, nem
n = j+j+m-2
call dlfkp(m, n, cp)
do i=1, nte
call dlftp (m, n, thet(i), cp, ped(i, j+ms2, ip))
end do
202 if (m>0) ped(1, j+ms2, ip) = 0.0
205 continue
!
else
!
do 207 j=1, nem
n = j+j+m-2
if (m>1.and.n>mxtr) then
do i=1, nte
u(i, j+ms2) = ped(i, j+ms2, ip)
end do
go to 207
end if
a1 = b(n-1)*a(n+m-3)/a(n+m-1)
b1 = a(n-m+1)/a(n+m-1)
if (n-m<=1) then
do i=1, nte
u(i, j+ms2) = a1*ped(i, j+ms2-1, ip) &
                   - b1*ped(i, j+ms2, ip)    
end do
else
c1 = b(n-1)*a(n-m-1)/a(n+m-1)
do i=1, nte
u(i, j+ms2) = a1*ped(i, j+ms2-1, ip) &
   - b1*ped(i, j+ms2, ip) + c1*u(i, j+ms2-1)    
end do
end if
207 continue
do j=1, nem
do i=1, nte
ped(i, j+ms2, ip) = u(i, j+ms2)
end do
end do
end if
!
if (ms2<=0.or.ms2>=nte) go to 200
!
! initialize array with random numbers using 
! Fortran90 intrinsics RANDOM_{SEED, NUMBER}
!
! old code commented out
!     do i=1, nte
!     xx(i) = rand()
!     end do
!
! replacement code
!
call RANDOM_SEED()
call random_number(xx(1:nte))
it = 0
201 do i=1, nte
z(i) = 0.0
wx(i) = 0.0
do j=1, nte
wx(i) = wx(i)+we(i, j, ip)*xx(j)
end do
end do
do 220 j=1, nte
if (j==ms2) go to 220
call gs(nte, wx, ped(1, j, ip), z)
220 continue
!  
do i=1, nte
xx(i) = xx(i)-z(i)
end do
call normal(nte, xx, idp, we(1, 1, ip))
it = it+1
if (it<=2) go to 201
do i=1, nte
ped(i, ms2, ip) = xx(i)
end do
200 continue
!
!     reorder if mtrunc is less than nlat-1 
!         case of even functions
!
if (modn==0) then
nshe(1) = (nlat-mtrunc-1)/2
nshe(2) = (nlat-mtrunc-2)/2
else
nshe(1) = (nlat-mtrunc)/2
nshe(2) = (nlat-mtrunc-1)/2
end if
!
do 210 mp1=1, 2
do j=1, nte
js = j+nshe(mp1)
if (js>nte) js = js-nte
do i=1, nte
u(i, js) = ped(i, j, mp1)
end do
end do
do j=1, nte
do i=1, nte
ped(i, j, mp1) = u(i, j)
end do
end do
210 continue
!
call trunc(0, nte, idp, ped(1, 1, 1), nte, ipse(1, 1))
call trunc(0, nte, idp, ped(1, 1, 2), nte, ipse(1, 2))
!
!     compute the analysis matrices
!
do 250 ip=1, 2
do i=1, nte
lock = 0
do j=1, nte
sum = 0.0
do k=1, nte
sum = sum+ped(k, i, ip)*we(k, j, ip)
end do
pe(i, j, ip) = ped(i, j, ip)
ze(j, i, ip) =  sum
if (abs(sum)>eps .and. lock==0) then
lock = 1
jzse(i, ip) = j
end if
end do
end do
250 continue
!
!     compute weight matrices for odd functions
!
do 50 mp1=1, 2
m = mp1-1
mrank = nlat-m-m
nem = (mrank+1)/2
nom = mrank-nem
do j=1, nom
n = j+j+m-1
call dlfkp(m, n, cp)
do i=1, nte
call dlftp (m, n, thet(i), cp, pod(i, j, mp1))
end do
if (modn==1) pod(nte, j, mp1) = 0.0
end do
call dsvdc(pod(m+1, 1, mp1), idp, nom, nom, s, e, u, &
                        idp, v(1,1), idp, work, 10, info)
!
do j=1, nom
s(j) = 1.0/(s(j)*s(j)) 
end do
!
!     compute weight matrix as u  s sup -2 u transpose
!
do j=1, nte
do i=1, nte
wo(i, j, mp1) = 0.0
end do
end do
do i=1, nom
do j=1, nom
sum = 0.
do k=1, nom
sum = sum+s(k)*u(i, k)*u(j, k)
end do
wo(i+m, j+m, mp1) = sum
end do
end do
50 continue
wo(1, 1, 2) = 1.0  
if (modn==1) then
wo(nte, nte, 1) = 1.0  
wo(nte, nte, 2) = 1.0  
end if
!
!     compute n**2 basis (odd functions)
!
ip = 2
do 300 mp1=1, mxtr+1
ip = 3-ip
m = mp1-1
ms2 = mp1/2
nrank = ms2+ms2
mrank = nlat-nrank
nem = (mrank+1)/2
nom = mrank-nem
!
!     compute associated legendre functions
!
if (m<=1) then
do 305 j=1, nom
n = j+j+m-1
call dlfkp(m, n, cp)
do i=1, nte
call dlftp (m, n, thet(i), cp, pod(i, j+ms2, ip))
end do
302 if (modn==1) pod(nte, j+ms2, ip) = 0.0
if (m>0) pod(1, j+ms2, ip) = 0.0
305 continue
!
else
!
do 307 j=1, nom
n = j+j+m-1
if (m>1.and.n>mxtr) then
do i=1, nte
u(i, j+ms2) = pod(i, j+ms2, ip)
end do
go to 304
end if
a1 = b(n-1)*a(n+m-3)/a(n+m-1)
b1 = a(n-m+1)/a(n+m-1)
if (n-m<=1) then
do i=1, nte
u(i, j+ms2) = a1*pod(i, j+ms2-1, ip) &
                   - b1*pod(i, j+ms2, ip)    
end do
else
c1 = b(n-1)*a(n-m-1)/a(n+m-1)
do i=1, nte
u(i, j+ms2) = a1*pod(i, j+ms2-1, ip) &
   - b1*pod(i, j+ms2, ip) + c1*u(i, j+ms2-1)    
end do
end if
304 if (modn==1) u(nte, j+ms2) = 0.0
307 continue
do j=1, nom
do i=1, nte
pod(i, j+ms2, ip) = u(i, j+ms2)
end do
end do
end if
!
if (ms2<=0.or.ms2>=nto) go to 300
!
! initialize array with random numbers using 
! Fortran90 intrinsics RANDOM_{SEED, NUMBER}
!
! old code commented out
!
!     do i=1, nte
!     xx(i) = rand()
!     end do
! replacement code
!
call random_number(xx(1:nte))
if (modn==1) xx(nte) = 0.0
it = 0
306 do i=1, nte
z(i) = 0.
wx(i) = 0.
do j=1, nto
wx(i) = wx(i)+wo(i, j, ip)*xx(j)
end do
end do
do 330 j=1, nto
if (j==ms2) go to 330
call gs(nte, wx, pod(1, j, ip), z(1))
330 continue
!  
do i=1, nte
xx(i) = xx(i)-z(i)
end do
call normal(nte, xx, idp, wo(1, 1, ip))
it = it+1
if (it<=2) go to 306
do i=1, nte
pod(i, ms2, ip) = xx(i)
end do
if (modn==1) pod(nte, ms2, ip) = 0.0
300 continue
!
!     reorder if mtrunc is less than nlat-1
!        case of odd functions  
!
if (modn==0) then
nsho(1) = (nlat-mtrunc)/2
nsho(2) = (nlat-mtrunc-1)/2
else
nsho(1) = (nlat-mtrunc-1)/2
nsho(2) = (nlat-mtrunc-2)/2
end if
!
do 310 mp1=1, 2
do j=1, nto
js = j+nsho(mp1)
if (js>nto) js = js-nto
do i=1, nte
u(i, js) = pod(i, j, mp1)
end do
end do
do j=1, nto
do i=1, nte
pod(i, j, mp1) = u(i, j)
end do
end do
310 continue
!
call trunc(0, nte, idp, pod(1, 1, 1), nto, ipso(1, 1))
call trunc(0, nte, idp, pod(1, 1, 2), nto, ipso(1, 2))
!
!     compute the analysis matrices (odd functions)
!
do ip=1, 2
do i=1, nto
lock = 0
do j=1, nto
sum = 0.0  
do k=1, nte
sum = sum+pod(k, i, ip)*wo(k, j, ip)
end do
po(i, j, ip) = pod(i, j, ip)
zo(j, i, ip) = sum
if (abs(sum)>eps .and. lock==0) then
lock = 1
jzso(i, ip) = j
end if
end do
end do
end do
return
end subroutine shpei1
!
!
! ... file shpe.f
!
! ... files which must be loaded with shpe.f
!
!     sphcom.f, hrfft.f
!
!     the n**2 projection with complement, odd/even
!     factorization and zero truncation on an
!     equally spaced grid as defined in the JCP paper
!     "Generalized discrete spherical harmonic transforms" 
!     by Paul N. Swarztrauber and William F. Spotz
!     It is equivalent to a harmonic analysis followed
!     by a synthesis except faster and requires less memory.
!
!     subroutine shpe(nlat, nlon, isym, mtrunc, x, y, idxy, 
!    1        wshp, lwshp, iwshp, liwshp, work, lwork, ierror)
!
!     shpe projects the array x onto the set of functions represented
!     by a discrete set of spherical harmonics.
!
!     input parameters
!
!     nlat   the number of colatitudes on the full sphere including the
!            poles. for example, nlat = 37 for a five degree grid.
!            nlat determines the grid increment in colatitude as
!            pi/(nlat-1).  if nlat is odd the equator is located at
!            grid point i=(nlat+1)/2. if nlat is even the equator is
!            located half way between points i=nlat/2 and i=nlat/2+1.
!            nlat must be at least 3.
!
!     nlon   the number of distinct londitude points.  nlon determines
!            the grid increment in longitude as 2*pi/nlon. for example
!            nlon = 72 for a five degree grid. nlon must be greater
!            than or equal to 4. the efficiency of the computation is
!            improved when nlon is a product of small prime numbers.
!            nlon must beat least 4.
!
!     isym   currently not used.    
!
!     mtrunc the highest longitudinal wave number retained in the
!            projection. It must be less than or equal to
!            the minimum of nlat-1 and nlon/2. The first wave
!            number is zero. For example, if wave numbers 0 and
!            1 are desired then mtrunc = 1.

!            zero.
!
!     x      a two dimensional array that contains the the nlat
!            by nlon array x(i, j) defined at the colatitude point 
!            theta(i) = (i-1)*pi/(nlat-1) and longitude point phi(j) =
!            (j-1)*2*pi/nlon.
!
!     idxy   the first dimension of the arrays x and y as they
!            appear in the program that calls shpe. It must be
!            at least nlat. 
!
!     wshp   a single precision array that must be saved for
!            repeated use by subroutine shpe.        
!
!     lwshp  the dimension of the array wshp as it appears in the
!            program that calls shpei. It must be at least
!            2*(nlat+1)**2+nlon+log2(nlon)
!
!     iwshp  an integer array that must be saved for repeated
!            use by subroutine shpe.        
!
!
!     liwshp the dimension of the array iwshp as it appears in the
!            program that calls shpei. It must be at least
!            4*(nlat+1).
!
!     work   a single precision work array that does 
!            not have to be saved.
!
!     lwork  the dimension of the array work as it appears in the
!            program that calls shpe. It must be at least
!            max(nlat*nlon, 4*(nlat+1)).
!
!     **************************************************************
!
!     output parameters
!
!     y      an nlat by nlon single precision array that contains 
!            the projection of x onto the set of functions that
!            can be represented by the discrete set of spherical 
!            harmonics. The arrays x(i, j) and y(i, j) are located
!            at colatitude point theta(i) = (i-1)*pi/(nlat-1) and 
!            longitude point phi(j) = (j-1)*2*pi/nlon.
!
!     ierror = 0  no errors
!            = 1  error in the specification of nlat
!            = 2  error in the specification of nlon
!            = 3  error in the specification of isym
!            = 4  error in the specification of mtrunc
!            = 5  error in the specification of lwshp
!            = 6  error in the specification of liwshp
!            = 7  error in the specification of lwork
!
subroutine shpe(nlat, nlon, isym, mtrunc, x, y, idxy, &
        wshp, lwshp, iwshp, liwshp, work, lwork, ierror)
!
dimension wshp(*), iwshp(*), work(*), x(idxy, nlon), y(idxy, nlon)
!
ierror = 1
if (nlat < 3) return
ierror = 2
if (nlon < 4) return
!      ierror = 3
!      if (isym.lt.0 .or. isym.gt.2) return
ierror = 4
mmax = min(nlat-1, nlon/2)
if (mtrunc<0 .or. mtrunc>mmax) return
ierror = 5
log2n = log(real(nlon))/log(2.0)
lw1 = 2*(nlat+1)**2
if (lwshp<lw1+nlon+log2n) return
ierror = 6
if (liwshp<4*(nlat+1)) return
ierror = 7
mwrk = max(nlat*nlon, 4*(nlat+1))
if (lwork <mwrk) return
ierror = 0
!
do j=1, nlon
 do i=1, nlat
  y(i, j) = x(i, j)
 end do
end do
call hrfftf(nlat, nlon, y, idxy, wshp(lw1+1), work)
!
nte = (nlat+1)/2
nloc1 = 2*nte*nte
nloc2 = nlat+1
iw1 = 1
iw2 = iw1+nloc1
iw3 = iw2+nloc1
iw4 = iw3+nloc1
jw1 = 1
jw2 = jw1+nloc2
jw3 = jw2+nloc2
jw4 = jw3+nloc2
!
call shpe1(nlat, nlon, isym, mtrunc, y, y, idxy, ierror, &
 nte, wshp(iw1), wshp(iw2), wshp(iw3), wshp(iw4), iwshp(jw1), &
 iwshp(jw2), iwshp(jw3), iwshp(jw4), work(jw1), &
 work(jw2), work(jw3), work(jw4))
!
call hrfftb(nlat, nlon, y, idxy, wshp(lw1+1), work)
!
sn = 1.0/nlon
do j=1, nlon
 do i=1, nlat
  y(i, j) = sn*y(i, j)
 end do
end do
return
end subroutine shpe
subroutine shpe1(nlat, nlon, isym, mtrunc, sx, sy, idxy, ierror, &
 idp, pe, po, ze, zo, ipse, jzse, ipso, jzso, xe, xo, ye, yo)
!
dimension sx(idxy, nlon), sy(idxy, nlon), nshe(2), nsho(2), &
  pe(idp, idp, 2), po(idp, idp, 2), ze(idp, idp, 2), zo(idp, idp, 2), &
  ipse(idp, 2), jzse(idp, 2), ipso(idp, 2), jzso(idp, 2), &
  xe(idp, 2), xo(idp, 2), ye(idp, 2), yo(idp, 2)
!
ns2 = nlat/2
modn = nlat-ns2-ns2
nte = (nlat+1)/2
nto = nlat-nte
!
if (modn==0) then
nshe(1) = (nlat-mtrunc-1)/2
nshe(2) = (nlat-mtrunc-2)/2
nsho(1) = (nlat-mtrunc)/2
nsho(2) = (nlat-mtrunc-1)/2
else
nshe(1) = (nlat-mtrunc)/2
nshe(2) = (nlat-mtrunc-1)/2
nsho(1) = (nlat-mtrunc-1)/2
nsho(2) = (nlat-mtrunc-2)/2
end if
mxtr = min(nlat-1, nlon/2, mtrunc)
ip = 2 
do 100 mp1=1, mxtr+1
ip = 3-ip
if (mxtr==nlat-1.and.mp1<=2) then
do i=1, nlat
sy(i, mp1) = sx(i, mp1)
end do
if (mp1==2) then
sy(1, 2) = 0.
sy(nlat, 2) = 0.
end if
if (nlon>=3) then
sy(1, 3) = 0.
sy(nlat, 3) = 0.
do i=2, nlat-1
sy(i, 3) = sx(i, 3)
end do
end if
go to 100
end if
m = mp1-1
mpm = max(1, m+m)
ms2 = mp1/2
mrank = min(nlat-m, nlat-ms2-ms2)   
!      mrank = mxtr+1-ms2-ms2
nrank = nlat-mrank
nem = (mrank+1)/2-nshe(ip)
nom = mrank-(mrank+1)/2-nsho(ip)
nec = nte-nem
noc = nto-nom
!
do i=1, nte
xe(i, 1) = .5*(sx(i, mpm)+sx(nlat+1-i, mpm))
xo(i, 1) = .5*(sx(i, mpm)-sx(nlat+1-i, mpm))
end do
if (mpm<nlon) then
do i=1, nte
xe(i, 2) = .5*(sx(i, mpm+1)+sx(nlat+1-i, mpm+1))
xo(i, 2) = .5*(sx(i, mpm+1)-sx(nlat+1-i, mpm+1))
end do
end if
if (3*nec<2*nem.or.nem==0) then
call tmxmx(nte, nec, idp, pe(1, 1, ip), nte, idp, &
          ze(1, 1, ip), xe, ye, ipse(1, ip), jzse(1, ip))  
do i=1, nte
ye(i, 1) = xe(i, 1)-ye(i, 1)
end do
if (mpm<nlon.and.m/=0) then
do i=1, nte
ye(i, 2) = xe(i, 2)-ye(i, 2)
end do
end if
else
call tmxmx(nte, nem, idp, pe(1, nec+1, ip), nte, idp, &
ze(1, nec+1, ip), xe, ye, ipse(nec+1, ip), jzse(nec+1, ip))
end if
if (3*noc<2*nom.or.nom==0) then
call tmxmx(nto, noc, idp, po(1, 1, ip), nto, idp, &
          zo(1, 1, ip), xo, yo, ipso(1, ip), jzso(1, ip))
do i=1, nte
yo(i, 1) = xo(i, 1)-yo(i, 1)
end do
if (mpm<nlon.and.m/=0) then
do i=1, nte
yo(i, 2) = xo(i, 2)-yo(i, 2)
end do
end if
else
call tmxmx(nto, nom, idp, po(1, noc+1, ip), nto, idp, &
zo(1, noc+1, ip), xo, yo, ipso(noc+1, ip), jzso(noc+1, ip))  
end if
do i=1, nte
sy(i, mpm) = ye(i, 1)+yo(i, 1)
sy(nlat+1-i, mpm) = ye(i, 1)-yo(i, 1)
end do
if (mpm<nlon.and.m/=0) then
do i=1, nte
sy(i, mpm+1) = ye(i, 2)+yo(i, 2)
sy(nlat+1-i, mpm+1) = ye(i, 2)-yo(i, 2)
end do 
end if
100 continue
js = mxtr+mxtr+2
do j=js, nlon
do i=1, nlat
sy(i, j) = 0.
end do
end do
return
end subroutine shpe1
subroutine mxm(lr, lc, ld, a, mc, md, b, nd, c)
real a(ld, *), b(md, *), c(nd, *)
do i=1, lr
do j=1, mc
c(i, j) = 0.
do k=1, lc 
c(i, j) = c(i, j)+a(i, k)*b(k, j)
end do
end do
end do
return
end subroutine mxm
subroutine smxm(lr, lc, ld, a, mc, md, b, nd, c)
dimension a(ld, *), b(md, *), c(nd, *)
do i=1, lr
do j=1, mc
c(i, j) = 0.
do k=1, lc 
c(i, j) = c(i, j)+a(i, k)*b(k, j)
end do
end do
end do
return
end subroutine smxm
subroutine mxmx(lr, lc, ld, a, mc, md, b, x, y)
dimension a(ld, *), b(md, *), x(ld, 2), y(ld, 2)
do k=1, lr
y(k, 1) = 0.
y(k, 2) = 0.
end do
!
if (lc<=0) return
do i=1, lc
sum1 = 0.
sum2 = 0.
do j=1, mc
sum1 = sum1 + b(i, j)*x(j, 1)
sum2 = sum2 + b(i, j)*x(j, 2)
end do
do k=1, lr
y(k, 1) = y(k, 1)+sum1*a(k, i)
y(k, 2) = y(k, 2)+sum2*a(k, i)
end do
end do
return
end subroutine mxmx
subroutine dmxmx(lr, lc, ld, a, mc, md, b, x, y)
real a(ld, *), b(md, *), x(ld, 2), y(ld, 2), &
                 sum1, sum2
do k=1, lr
y(k, 1) = 0.
y(k, 2) = 0.
end do
!
if (lc<=0) return
do i=1, lc
sum1 = 0.
sum2 = 0.
do j=1, mc
sum1 = sum1 + b(i, j)*x(j, 1)
sum2 = sum2 + b(i, j)*x(j, 2)
end do
do k=1, lr
y(k, 1) = y(k, 1)+sum1*a(k, i)
y(k, 2) = y(k, 2)+sum2*a(k, i)
end do
end do
return
end subroutine dmxmx
subroutine tmxmx(lr, lc, ld, a, mc, md, b, x, y, is, js)
dimension a(ld, *), b(md, *), x(ld, 2), y(ld, 2), &
              is(*), js(*)
!
kmx = min(lr+1, ld)
do k=1, kmx
y(k, 1) = 0.
y(k, 2) = 0.
end do
if (lc<=0) return
!
do i=1, lc
sum1 = 0.
sum2 = 0.
do j=js(i), mc
sum1 = sum1 + b(j, i)*x(j, 1)
sum2 = sum2 + b(j, i)*x(j, 2)
end do
do k=is(i), lr
y(k, 1) = y(k, 1)+sum1*a(k, i)
y(k, 2) = y(k, 2)+sum2*a(k, i)
end do
end do
return
end subroutine tmxmx
subroutine trunc(irc, n, idp, a, nrc, ijs)
real a, eps
parameter (eps=epsilon(1.0))
dimension a(idp, *), ijs(n)
!
!     irc = 0 for columns , or irc = 1 for rows
!
if (irc/=0) go to 30
do 20 j=1, nrc
do i=1, n
ijs(j) = i
if (abs(a(i, j)) > eps) go to 20
end do
20 continue
return
30 do 50 i=1, nrc
do j=1, n
ijs(i) = j
if (abs(a(i, j)) > eps) go to 50
end do
50 continue
return
end subroutine trunc
subroutine gs(n, x, y, z)
dimension x(n), y(n), z(n)
real x, y, z, sum
!
!     accumulate innerproducts of x with respect to y.
!
sum = 0.
do i=1, n
sum = sum+x(i)*y(i)
end do
do i=1, n
z(i) = z(i)+sum*y(i)
end do
return
end subroutine gs
subroutine normal(n, x, id, q)
dimension x(n), q(id, n)
real x, q, sum, sqs
!
!     normalize x
!
sqs = 0.
do i=1, n
sum = 0.
do j=1, n
sum = sum+q(i, j)*x(j)
end do
sqs = sqs+sum*x(i)
end do
!
sqs = sqrt(sqs)
do i=1, n
x(i) = x(i)/sqs
end do
return
end subroutine normal
subroutine coe(moe, n, x, dmax)
real x(n), dmax
nh = (n+1)/2
dmax = 0.
if (moe/=0) go to 1
do i=1, nh
dmax = max(dmax, abs(x(i)-x(n-i+1)))
x(i) = .5*(x(i)+x(n-i+1))
x(n-i+1) = x(i)
end do
return
1 do i=1, nh
dmax = max(dmax, abs(x(i)+x(n-i+1)))
x(i) = .5*(x(i)-x(n-i+1))
x(n-i+1) = -x(i)
end do
if (mod(n, 2)/=0) x(nh) = 0.
return
end subroutine coe
!     subroutine dlfkp(m, n, cp)
!
!     subroutine dlfkp computes the coefficients in the trigonometric
!     expansion of the normalized associated legendre functions:
!
!     pbar(m, n, theta) = sqrt((2*n+1)*factorial(n-m)/(2*factorial(n+m)))
!                        *sin(theta)**m/(2**n*factorial(n)) times the
!                        (n+m)th derivative of (x**2-1)**n with respect
!                        to x=cos(theta)
!
!     where theta is colatitude.
!
!     subroutine dlfkp computes the coefficients cp(k) in the
!     following trigonometric expansion of pbar(m, n, theta).
!
!            1) for n even and m even, pbar(m, n, theta) =
!               .5*cp(1) plus the sum from k=1 to k=n/2
!               of cp(k)*cos(2*k*th)
!
!            2) for n even and m odd, pbar(m, n, theta) =
!               the sum from k=1 to k=n/2 of
!               cp(k)*sin(2*k*th)
!
!            3) for n odd and m even, pbar(m, n, theta) =
!               the sum from k=1 to k=(n+1)/2 of
!               cp(k)*cos((2*k-1)*th)
!
!            4) for n odd and m odd,  pbar(m, n, theta) =
!               the sum from k=1 to k=(n+1)/2 of
!               cp(k)*sin((2*k-1)*th)
!
!     input parameters
!
!     m      is the order of pbar(n, m, theta). m can be any integer
!            however pbar(n, m, theta) = 0  if abs(m) is greater than
!            n and pbar(n, m, theta) = (-1)**m*pbar(n, -m, theta) for
!            negative m.
!
!     n      nonnegative integer specifying the degree of
!            pbar(n, m, theta)
!
!     output parameters
!
!     cp     a real array that contains the fourier
!            coefficients for pbar(m, n, theta). the length of the
!            array depends on the parity of m and n
!
!                  parity            length of cp
!
!               n even m even           n/2+1
!               n even m odd             n/2
!               n odd  m even          (n+1)/2
!               n odd  m odd           (n+1)/2
!
!
! ****************************************************************
subroutine dlfkp (m, n, cp)
!
real cp, fnum, fden, fnmh, a1, b1, c1, cp2, fnnp1, fnmsq, fk, &
       t1, t2, pm1, sc10, sc20, sc40
dimension       cp(1)
parameter (sc10=1024.0)
parameter (sc20=sc10*sc10)
parameter (sc40=sc20*sc20)
!
cp(1) = 0.
ma = abs(m)
if (ma > n) return
if (n-1< 0) then
    goto 2
else if (n-1 == 0) then 
    goto 3
else 
    goto 5
end if
2 cp(1) = sqrt(2.0)
return
3 if (ma /= 0) go to 4
cp(1) = sqrt(1.5)
return
4 cp(1) = sqrt(.75)
if (m == -1) cp(1) = -cp(1)
return
5 if (mod(n+ma, 2) /= 0) go to 10
nmms2 = (n-ma)/2
fnum = n+ma+1
fnmh = n-ma+1
pm1 = 1.0
go to 15
10 nmms2 = (n-ma-1)/2
fnum = n+ma+2
fnmh = n-ma+2
pm1 = -1.0
!      t1 = 1.
!      t1 = 2.0**(n-1)
!      t1 = 1.0/t1
15 t1 = 1.0/sc20
nex = 20
fden = 2.0
if (nmms2 < 1) go to 20
do 18 i=1, nmms2
t1 = fnum*t1/fden
if (t1 > sc20) then
t1 = t1/sc40
nex = nex+40
end if
fnum = fnum+2.
fden = fden+2.
18 continue
20 t1 = t1/2.0**(n-1-nex)
if (mod(ma/2, 2) /= 0) t1 = -t1
t2 = 1. 
if (ma == 0) go to 26
do 25 i=1, ma
t2 = fnmh*t2/(fnmh+pm1)
fnmh = fnmh+2.
25 continue
26 cp2 = t1*sqrt((n+.5)*t2)
fnnp1 = n*(n+1)
fnmsq = fnnp1-2.0*ma*ma
l = (n+1)/2
if (mod(n, 2) == 0 .and. mod(ma, 2) == 0) l = l+1
cp(l) = cp2
if (m >= 0) go to 29
if (mod(ma, 2) /= 0) cp(l) = -cp(l)
29 if (l <= 1) return
fk = n
a1 = (fk-2.)*(fk-1.)-fnnp1
b1 = 2.*(fk*fk-fnmsq)
cp(l-1) = b1*cp(l)/a1
30 l = l-1
if (l <= 1) return
fk = fk-2.
a1 = (fk-2.)*(fk-1.)-fnnp1
b1 = -2.*(fk*fk-fnmsq)
c1 = (fk+1.)*(fk+2.)-fnnp1
cp(l-1) = -(b1*cp(l)+c1*cp(l+1))/a1
go to 30
end subroutine dlfkp
subroutine dlftp (m, n, theta, cp, pb)
dimension cp(1)
real cp, pb, theta, cdt, sdt, cth, sth, chh
cdt = cos(2.0*theta)
sdt = sin(2.0*theta)
nmod=mod(n, 2)
mmod=mod(abs(m), 2)
if (nmod< 0) then
    goto 1
else if (nmod == 0) then 
    goto 1
else 
    goto 2
end if
1 if (mmod< 0) then
    goto 3
else if (mmod == 0) then 
    goto 3
else 
    goto 4
end if
!
!     n even, m even
!
3 kdo=n/2
pb = .5*cp(1)
if (n == 0) return
cth = cdt
sth = sdt
do 170 k=1, kdo
!     pb = pb+cp(k+1)*cos(2*k*theta)
pb = pb+cp(k+1)*cth
chh = cdt*cth-sdt*sth
sth = sdt*cth+cdt*sth
cth = chh
170 continue
return
!
!     n even, m odd
!
4 kdo = n/2
pb = 0.
cth = cdt
sth = sdt
do 180 k=1, kdo
!     pb = pb+cp(k)*sin(2*k*theta)
pb = pb+cp(k)*sth
chh = cdt*cth-sdt*sth
sth = sdt*cth+cdt*sth
cth = chh
180 continue
return
2 if (mmod< 0) then
    goto 13
else if (mmod == 0) then 
    goto 13
else 
    goto 14
end if
!
!     n odd, m even
!
13 kdo = (n+1)/2
pb = 0.
cth = cos(theta)
sth = sin(theta)
do 190 k=1, kdo
!     pb = pb+cp(k)*cos((2*k-1)*theta)
pb = pb+cp(k)*cth
chh = cdt*cth-sdt*sth
sth = sdt*cth+cdt*sth
cth = chh
190 continue
return
!
!     n odd, m odd
!
14 kdo = (n+1)/2
pb = 0.
cth = cos(theta)
sth = sin(theta)
do 200 k=1, kdo
!     pb = pb+cp(k)*sin((2*k-1)*theta)
pb = pb+cp(k)*sth
chh = cdt*cth-sdt*sth
sth = sdt*cth+cdt*sth
cth = chh
200 continue
return
end subroutine dlftp
!
subroutine dsvdc(x, ldx, n, p, s, e, u, ldu, v, ldv, work, job, info)
integer ldx, n, p, ldu, ldv, job, info
real x(ldx, 1), s(1), e(1), u(ldu, 1), v(ldv, 1), work(1)
!
!
!     dsvdc is a subroutine to reduce a real nxp matrix x
!     by orthogonal transformations u and v to diagonal form.  the
!     diagonal elements s(i) are the singular values of x.  the
!     columns of u are the corresponding left singular vectors, 
!     and the columns of v the right singular vectors.
!
!     on entry
!
!         x         real(ldx, p), where ldx.ge.n.
!                   x contains the matrix whose singular value
!                   decomposition is to be computed.  x is
!                   destroyed by dsvdc.
!
!         ldx       integer.
!                   ldx is the leading dimension of the array x.
!
!         n         integer.
!                   n is the number of rows of the matrix x.
!
!         p         integer.
!                   p is the number of columns of the matrix x.
!
!         ldu       integer.
!                   ldu is the leading dimension of the array u.
!                   (see below).
!
!         ldv       integer.
!                   ldv is the leading dimension of the array v.
!                   (see below).
!
!         work      real(n).
!                   work is a scratch array.
!
!         job       integer.
!                   job controls the computation of the singular
!                   vectors.  it has the decimal expansion ab
!                   with the following meaning
!
!                        a.eq.0    do not compute the left singular
!                                  vectors.
!                        a.eq.1    return the n left singular vectors
!                                  in u.
!                        a.ge.2    return the first min(n, p) singular
!                                  vectors in u.
!                        b.eq.0    do not compute the right singular
!                                  vectors.
!                        b.eq.1    return the right singular vectors
!                                  in v.
!
!     on return
!
!         s         real(mm), where mm=min(n+1, p).
!                   the first min(n, p) entries of s contain the
!                   singular values of x arranged in descending
!                   order of magnitude.
!
!         e         real(p), 
!                   e ordinarily contains zeros.  however see the
!                   discussion of info for exceptions.
!
!         u         real(ldu, k), where ldu.ge.n.  if
!                                   joba.eq.1 then k.eq.n, if joba.ge.2
!                                   then k.eq.min(n, p).
!                   u contains the matrix of left singular vectors.
!                   u is not referenced if joba.eq.0.  if n.le.p
!                   or if joba.eq.2, then u may be identified with x
!                   in the subroutine call.
!
!         v         real(ldv, p), where ldv.ge.p.
!                   v contains the matrix of right singular vectors.
!                   v is not referenced if job.eq.0.  if p.le.n, 
!                   then v may be identified with x in the
!                   subroutine call.
!
!         info      integer.
!                   the singular values (and their corresponding
!                   singular vectors) s(info+1), s(info+2), ..., s(m)
!                   are correct (here m=min(n, p)).  thus if
!                   info.eq.0, all the singular values and their
!                   vectors are correct.  in any event, the matrix
!                   b = trans(u)*x*v is the bidiagonal matrix
!                   with the elements of s on its diagonal and the
!                   elements of e on its super-diagonal (trans(u)
!                   is the transpose of u).  thus the singular
!                   values of x and b are the same.
!
!     linpack. this version dated 08/14/78 .
!              correction made to shift 2/84.
!     g.w. stewart, university of maryland, argonne national lab.
!
!     dsvdc uses the following functions and subprograms.
!
!     external drot
!     blas daxpy, ddot, dscal, dswap, dnrm2, drotg
!     fortran dabs, dmax1, max0, min0, mod, dsqrt
!
!     internal variables
!
integer i, iter, j, jobu, k, kase, kk, l, ll, lls, lm1, lp1, ls, lu, m, maxit, &
        mm, mm1, mp1, nct, nctp1, ncu, nrt, nrtp1
!      real ddot, t, r
real ddot, t
real b, c, cs, el, emm1, f, g, dnrm2, scale, shift, sl, sm, sn, &
                 smm1, t1, test, ztest
logical wantu, wantv
!
!
!     set the maximum number of iterations.
!
maxit = 30
!
!     determine what is to be computed.
!
wantu = .false.
wantv = .false.
jobu = mod(job, 100)/10
ncu = n
if (jobu > 1) ncu = min(n, p)
if (jobu /= 0) wantu = .true.
if (mod(job, 10) /= 0) wantv = .true.
!
!     reduce x to bidiagonal form, storing the diagonal elements
!     in s and the super-diagonal elements in e.
!
info = 0
nct = min(n-1, p)
nrt = max(0, min(p-2, n))
lu = max(nct, nrt)
if (lu < 1) go to 170
do 160 l = 1, lu
   lp1 = l + 1
   if (l > nct) go to 20
!
!           compute the transformation for the l-th column and
!           place the l-th diagonal in s(l).
!
      s(l) = dnrm2(n-l+1, x(l, l), 1)
      if (s(l) == 0.0) go to 10
         if (x(l, l) /= 0.0) s(l) = dsign(s(l), x(l, l))
         call dscal(n-l+1, 1.0/s(l), x(l, l), 1)
         x(l, l) = 1.0 + x(l, l)
10       continue
      s(l) = -s(l)
20    continue
   if (p < lp1) go to 50
   do 40 j = lp1, p
      if (l > nct) go to 30
      if (s(l) == 0.0) go to 30
!
!              apply the transformation.
!
         t = -ddot(n-l+1, x(l, l), 1, x(l, j), 1)/x(l, l)
         call daxpy(n-l+1, t, x(l, l), 1, x(l, j), 1)
30       continue
!
!           place the l-th row of x into  e for the
!           subsequent calculation of the row transformation.
!
      e(j) = x(l, j)
40    continue
50    continue
   if (.not.wantu .or. l > nct) go to 70
!
!           place the transformation in u for subsequent back
!           multiplication.
!
      do 60 i = l, n
         u(i, l) = x(i, l)
60       continue
70    continue
   if (l > nrt) go to 150
!
!           compute the l-th row transformation and place the
!           l-th super-diagonal in e(l).
!
      e(l) = dnrm2(p-l, e(lp1), 1)
      if (e(l) == 0.0) go to 80
         if (e(lp1) /= 0.0) e(l) = dsign(e(l), e(lp1))
         call dscal(p-l, 1.0/e(l), e(lp1), 1)
         e(lp1) = 1.0 + e(lp1)
80       continue
      e(l) = -e(l)
      if (lp1 > n .or. e(l) == 0.0) go to 120
!
!              apply the transformation.
!
         do 90 i = lp1, n
            work(i) = 0.0
90          continue
         do 100 j = lp1, p
            call daxpy(n-l, e(j), x(lp1, j), 1, work(lp1), 1)
100          continue
         do 110 j = lp1, p
            call daxpy(n-l, -e(j)/e(lp1), work(lp1), 1, x(lp1, j), 1)
110          continue
120       continue
      if (.not.wantv) go to 140
!
!              place the transformation in v for subsequent
!              back multiplication.
!
         do 130 i = lp1, p
            v(i, l) = e(i)
130          continue
140       continue
150    continue
160 continue
170 continue
!
!     set up the final bidiagonal matrix or order m.
!
m = min(p, n+1)
nctp1 = nct + 1
nrtp1 = nrt + 1
if (nct < p) s(nctp1) = x(nctp1, nctp1)
if (n < m) s(m) = 0.0
if (nrtp1 < m) e(nrtp1) = x(nrtp1, m)
e(m) = 0.0
!
!     if required, generate u.
!
if (.not.wantu) go to 300
   if (ncu < nctp1) go to 200
   do 190 j = nctp1, ncu
      do 180 i = 1, n
         u(i, j) = 0.0
180       continue
      u(j, j) = 1.0
190    continue
200    continue
   if (nct < 1) go to 290
   do 280 ll = 1, nct
      l = nct - ll + 1
      if (s(l) == 0.0) go to 250
         lp1 = l + 1
         if (ncu < lp1) go to 220
         do 210 j = lp1, ncu
            t = -ddot(n-l+1, u(l, l), 1, u(l, j), 1)/u(l, l)
            call daxpy(n-l+1, t, u(l, l), 1, u(l, j), 1)
210          continue
220          continue
         call dscal(n-l+1, -1.0, u(l, l), 1)
         u(l, l) = 1.0 + u(l, l)
         lm1 = l - 1
         if (lm1 < 1) go to 240
         do 230 i = 1, lm1
            u(i, l) = 0.0
230          continue
240          continue
      go to 270
250       continue
         do 260 i = 1, n
            u(i, l) = 0.0
260          continue
         u(l, l) = 1.0
270       continue
280    continue
290    continue
300 continue
!
!     if it is required, generate v.
!
if (.not.wantv) go to 350
   do 340 ll = 1, p
      l = p - ll + 1
      lp1 = l + 1
      if (l > nrt) go to 320
      if (e(l) == 0.0) go to 320
         do 310 j = lp1, p
            t = -ddot(p-l, v(lp1, l), 1, v(lp1, j), 1)/v(lp1, l)
            call daxpy(p-l, t, v(lp1, l), 1, v(lp1, j), 1)
310          continue
320       continue
      do 330 i = 1, p
         v(i, l) = 0.0
330       continue
      v(l, l) = 1.0
340    continue
350 continue
!
!     main iteration loop for the singular values.
!
mm = m
iter = 0
360 continue
!
!        quit if all the singular values have been found.
!
!     ...exit
   if (m == 0) go to 620
!
!        if too many iterations have been performed, set
!        flag and return.
!
   if (iter < maxit) go to 370
      info = m
!     ......exit
      go to 620
370    continue
!
!        this section of the program inspects for
!        negligible elements in the s and e arrays.  on
!        completion the variables kase and l are set as follows.
!
!           kase = 1     if s(m) and e(l-1) are negligible and l.lt.m
!           kase = 2     if s(l) is negligible and l.lt.m
!           kase = 3     if e(l-1) is negligible, l.lt.m, and
!                        s(l), ..., s(m) are not negligible (qr step).
!           kase = 4     if e(m-1) is negligible (convergence).
!
   do 390 ll = 1, m
      l = m - ll
!        ...exit
      if (l == 0) go to 400
      test = abs(s(l)) + abs(s(l+1))
      ztest = test + abs(e(l))
      if (ztest /= test) go to 380
         e(l) = 0.0
!        ......exit
         go to 400
380       continue
390    continue
400    continue
   if (l /= m - 1) go to 410
      kase = 4
   go to 480
410    continue
      lp1 = l + 1
      mp1 = m + 1
      do 430 lls = lp1, mp1
         ls = m - lls + lp1
!           ...exit
         if (ls == l) go to 440
         test = 0.0
         if (ls /= m) test = test + abs(e(ls))
         if (ls /= l + 1) test = test + abs(e(ls-1))
         ztest = test + abs(s(ls))
         if (ztest /= test) go to 420
            s(ls) = 0.0
!           ......exit
            go to 440
420          continue
430       continue
440       continue
      if (ls /= l) go to 450
         kase = 3
      go to 470
450       continue
      if (ls /= m) go to 460
         kase = 1
      go to 470
460       continue
         kase = 2
         l = ls
470       continue
480    continue
   l = l + 1
!
!        perform the task indicated by kase.
!
   go to (490, 520, 540, 570), kase
!
!        deflate negligible s(m).
!
490    continue
      mm1 = m - 1
      f = e(m-1)
      e(m-1) = 0.0
      do 510 kk = l, mm1
         k = mm1 - kk + l
         t1 = s(k)
         call drotg(t1, f, cs, sn)
         s(k) = t1
         if (k == l) go to 500
            f = -sn*e(k-1)
            e(k-1) = cs*e(k-1)
500          continue
         if (wantv) call drot(p, v(1, k), 1, v(1, m), 1, cs, sn)
510       continue
   go to 610
!
!        split at negligible s(l).
!
520    continue
      f = e(l-1)
      e(l-1) = 0.0
      do 530 k = l, m
         t1 = s(k)
         call drotg(t1, f, cs, sn)
         s(k) = t1
         f = -sn*e(k)
         e(k) = cs*e(k)
         if (wantu) call drot(n, u(1, k), 1, u(1, l-1), 1, cs, sn)
530       continue
   go to 610
!
!        perform one qr step.
!
540    continue
!
!           calculate the shift.
!
      scale = dmax1(abs(s(m)), abs(s(m-1)), abs(e(m-1)), &
                    abs(s(l)), abs(e(l)))
      sm = s(m)/scale
      smm1 = s(m-1)/scale
      emm1 = e(m-1)/scale
      sl = s(l)/scale
      el = e(l)/scale
      b = ((smm1 + sm)*(smm1 - sm) + emm1**2)/2.0
      c = (sm*emm1)**2
      shift = 0.0
      if (b == 0.0 .and. c == 0.0) go to 550
         shift = sqrt(b**2+c)
         if (b < 0.0) shift = -shift
         shift = c/(b + shift)
550       continue
      f = (sl + sm)*(sl - sm) + shift
      g = sl*el
!
!           chase zeros.
!
      mm1 = m - 1
      do 560 k = l, mm1
         call drotg(f, g, cs, sn)
         if (k /= l) e(k-1) = f
         f = cs*s(k) + sn*e(k)
         e(k) = cs*e(k) - sn*s(k)
         g = sn*s(k+1)
         s(k+1) = cs*s(k+1)
         if (wantv) call drot(p, v(1, k), 1, v(1, k+1), 1, cs, sn)
         call drotg(f, g, cs, sn)
         s(k) = f
         f = cs*e(k) + sn*s(k+1)
         s(k+1) = -sn*e(k) + cs*s(k+1)
         g = sn*e(k+1)
         e(k+1) = cs*e(k+1)
         if (wantu .and. k < n) &
            call drot(n, u(1, k), 1, u(1, k+1), 1, cs, sn)
560       continue
      e(m-1) = f
      iter = iter + 1
   go to 610
!
!        convergence.
!
570    continue
!
!           make the singular value  positive.
!
      if (s(l) >= 0.0) go to 580
         s(l) = -s(l)
         if (wantv) call dscal(p, -1.0, v(1, l), 1)
580       continue
!
!           order the singular value.
!
590       if (l == mm) go to 600
!           ...exit
         if (s(l) >= s(l+1)) go to 600
         t = s(l)
         s(l) = s(l+1)
         s(l+1) = t
         if (wantv .and. l < p) &
            call dswap(p, v(1, l), 1, v(1, l+1), 1)
         if (wantu .and. l < n) &
            call dswap(n, u(1, l), 1, u(1, l+1), 1)
         l = l + 1
      go to 590
600       continue
      iter = 0
      m = m - 1
610    continue
go to 360
620 continue
return
end subroutine dsvdc
subroutine daxpy(n, da, dx, incx, dy, incy)
!
!     constant times a vector plus a vector.
!     uses unrolled loops for increments equal to one.
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
real dx(*), dy(*), da
integer i, incx, incy, ix, iy, m, mp1, n
!
if (n<=0)return
if (da == 0.0) return
if (incx==1.and.incy==1)go to 20
!
!        code for unequal increments or equal increments
!          not equal to 1
!
ix = 1
iy = 1
if (incx<0)ix = (-n+1)*incx + 1
if (incy<0)iy = (-n+1)*incy + 1
do 10 i = 1, n
  dy(iy) = dy(iy) + da*dx(ix)
  ix = ix + incx
  iy = iy + incy
10 continue
return
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
20 m = mod(n, 4)
if ( m == 0 ) go to 40
do 30 i = 1, m
  dy(i) = dy(i) + da*dx(i)
30 continue
if ( n < 4 ) return
40 mp1 = m + 1
do 50 i = mp1, n, 4
  dy(i) = dy(i) + da*dx(i)
  dy(i + 1) = dy(i + 1) + da*dx(i + 1)
  dy(i + 2) = dy(i + 2) + da*dx(i + 2)
  dy(i + 3) = dy(i + 3) + da*dx(i + 3)
50 continue

end subroutine daxpy


    pure function ddot(n, dx, incx, dy, incy) result (return_value)
        implicit none
        !
        !     forms the dot product of two vectors.
        !     uses unrolled loops for increments equal to one.
        !     jack dongarra, linpack, 3/11/78.
        !     modified 12/3/93, array(1) declarations changed to array(*)
        !
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        integer, intent (in) :: n
        real,    intent (in) :: dx(*)
        integer, intent (in) :: incx
        real,    intent (in) :: dy(*)
        integer, intent (in) :: incy
        real                 :: return_value
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        real    :: dtemp
        integer :: i, ix, iy, m, mp1
        !----------------------------------------------------------------------

        ! Initialize
        return_value = 0.0
        dtemp = 0.0

        if (n <= 0) then
            return
        end if

        if (incx == 1 .and. incy == 1) then
            go to 20
        end if
        !
        !        code for unequal increments or equal increments
        !          not equal to 1
        !
        ix = 1
        iy = 1
        if (incx < 0) then
            ix = (-n+1)*incx + 1
        end if

        if (incy < 0) then
            iy = (-n+1)*incy + 1
        end if

        do i = 1, n
            dtemp = dtemp + dx(ix)*dy(iy)
            ix = ix + incx
            iy = iy + incy
        end do

        return_value = dtemp
        return
        !
        !        code for both increments equal to 1
        !
        !
        !        clean-up loop
        !
20      m = mod(n, 5)
        if (m == 0) then
            go to 40
        end if

        do i = 1, m
            dtemp = dtemp + dx(i)*dy(i)
        end do

        if ( n < 5 ) then
            go to 60
        end if


40      mp1 = m + 1
        do i = mp1, n, 5
            dtemp = dtemp &
                + dx(i)*dy(i) &
                + dx(i + 1)*dy(i + 1) &
                + dx(i + 2)*dy(i + 2) &
                + dx(i + 3)*dy(i + 3)&
                + dx(i + 4)*dy(i + 4)
        end do
60      return_value = dtemp

    end function ddot


    pure function dnrm2( n, x, incx ) result (return_value)
        implicit none
        !
        ! Purpose:
        !
        !  dnrm2 returns the euclidean norm of a vector via the function
        !  name, so that
        !
        !     dnrm2 := sqrt( x'*x )
        !
        !
        !
        !  -- This version written on 25-October-1982.
        !     Modified on 14-October-1993 to inline the call to DLASSQ.
        !     Sven Hammarling, Nag Ltd.
        !
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        integer, intent (in) :: n
        real,    intent (in) :: x(*)
        integer, intent (in) :: incx
        real                 :: return_value
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer         :: ix
        real            :: absxi, norm, scale_rename, ssq
        real, parameter :: ZERO = nearest(1.0,1.0)-nearest(1.0,-1.0)
        real, parameter :: ONE = 1.0
        !----------------------------------------------------------------------

        !
        !==> Executable Statements
        !
        if (n < 1 .or. incx < 1 ) then
            norm = ZERO
        else if (n == 1) then
            norm = abs(x(1))
        else
            scale_rename = ZERO
            ssq = ONE
            !        The following loop is equivalent to this call to the LAPACK
            !        auxiliary routine:
            !        call dlassq( n, x, incx, scale, ssq )
            !
            do ix = 1, 1 + (n - 1) * incx, incx
                if ( x(ix) /= ZERO ) then
                    absxi = abs(x(ix))
                    if ( scale_rename < absxi ) then
                        ssq = ONE + ssq * (scale_rename/absxi)**2
                        scale_rename = absxi
                    else
                        ssq = ssq + (absxi/scale_rename)**2
                    end if
                end if
            end do
            norm  = scale_rename * sqrt( ssq )
        end if

        !
        !==> Return norm
        !
        return_value = norm

    end function dnrm2



subroutine  drot (n, dx, incx, dy, incy, c, s)
!
!     applies a plane rotation.
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
real dx(*), dy(*), dtemp, c, s
integer i, incx, incy, ix, iy, n
!
if (n<=0)return
if (incx==1.and.incy==1)go to 20
!
!       code for unequal increments or equal increments not equal
!         to 1
!
ix = 1
iy = 1
if (incx<0)ix = (-n+1)*incx + 1
if (incy<0)iy = (-n+1)*incy + 1
do 10 i = 1, n
  dtemp = c*dx(ix) + s*dy(iy)
  dy(iy) = c*dy(iy) - s*dx(ix)
  dx(ix) = dtemp
  ix = ix + incx
  iy = iy + incy
10 continue
return
!
!       code for both increments equal to 1
!
20 do 30 i = 1, n
  dtemp = c*dx(i) + s*dy(i)
  dy(i) = c*dy(i) - s*dx(i)
  dx(i) = dtemp
30 continue
return
end subroutine  drot
subroutine drotg(da, db, c, s)
!
!     construct givens plane rotation.
!     jack dongarra, linpack, 3/11/78.
!
real da, db, c, s, roe, scale, r, z
!
roe = db
if ( abs(da) > abs(db) ) roe = da
scale = abs(da) + abs(db)
if ( scale /= 0.0 ) go to 10
   c = 1.0
   s = 0.0
   r = 0.0
   z = 0.0
   go to 20
10 r = scale*sqrt((da/scale)**2 + (db/scale)**2)
r = dsign(1.0, roe)*r
c = da/r
s = db/r
z = 1.0
if ( abs(da) > abs(db) ) z = s
if ( abs(db) >= abs(da) .and. c /= 0.0 ) z = 1.0/c
20 da = r
db = z
return
end subroutine drotg
subroutine  dscal(n, da, dx, incx)
!
!     scales a vector by a constant.
!     uses unrolled loops for increment equal to one.
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
real da, dx(*)
integer i, incx, m, mp1, n, nincx
!
if ( n<=0 .or. incx<=0 )return
if (incx==1)go to 20
!
!        code for increment not equal to 1
!
nincx = n*incx
do 10 i = 1, nincx, incx
  dx(i) = da*dx(i)
10 continue
return
!
!        code for increment equal to 1
!
!
!        clean-up loop
!
20 m = mod(n, 5)
if ( m == 0 ) go to 40
do 30 i = 1, m
  dx(i) = da*dx(i)
30 continue
if ( n < 5 ) return
40 mp1 = m + 1
do 50 i = mp1, n, 5
  dx(i) = da*dx(i)
  dx(i + 1) = da*dx(i + 1)
  dx(i + 2) = da*dx(i + 2)
  dx(i + 3) = da*dx(i + 3)
  dx(i + 4) = da*dx(i + 4)
50 continue
return
end subroutine  dscal
subroutine  dswap (n, dx, incx, dy, incy)
!
!     interchanges two vectors.
!     uses unrolled loops for increments equal one.
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
real dx(*), dy(*), dtemp
integer i, incx, incy, ix, iy, m, mp1, n
!
if (n<=0)return
if (incx==1.and.incy==1)go to 20
!
!       code for unequal increments or equal increments not equal
!         to 1
!
ix = 1
iy = 1
if (incx<0)ix = (-n+1)*incx + 1
if (incy<0)iy = (-n+1)*incy + 1
do 10 i = 1, n
  dtemp = dx(ix)
  dx(ix) = dy(iy)
  dy(iy) = dtemp
  ix = ix + incx
  iy = iy + incy
10 continue
return
!
!       code for both increments equal to 1
!
!
!       clean-up loop
!
20 m = mod(n, 3)
if ( m == 0 ) go to 40
do 30 i = 1, m
  dtemp = dx(i)
  dx(i) = dy(i)
  dy(i) = dtemp
30 continue
if ( n < 3 ) return
40 mp1 = m + 1
do 50 i = mp1, n, 3
  dtemp = dx(i)
  dx(i) = dy(i)
  dy(i) = dtemp
  dtemp = dx(i + 1)
  dx(i + 1) = dy(i + 1)
  dy(i + 1) = dtemp
  dtemp = dx(i + 2)
  dx(i + 2) = dy(i + 2)
  dy(i + 2) = dtemp
50 continue
return
end subroutine  dswap
