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
!     *       A Package of Fortran Subroutines and Programs           *
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
!                           August 2003
!
! ... in file shpg.f
!
!     this file contains code and documentation for subroutines
!     shpgi and shpg.
!
! ... files which must be loaded with shpg.f
!
!     type_HFFTpack.f
!
!     shpgi initializes the arrays wshp and iwshp for subsequent 
!     use in subroutine shpg, which performs the harmonic projection 
!     which is equivalent to a harmonic analysis followed by 
!     harmonic synthesis but faster and with less memory.
!     (see description of subroutine shpg below).
!
!     subroutine shpgi(nlat, nlon, isym, mtrunc, wshp, lwshp, iwshp, 
!    1 liwshp, work, lwork, ierror)
!
!     shpgi initializes arrays wshp and iwshp for repeated use
!     by subroutine shpg ....
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
!            nlon must be at least 4.
!
!     isym   currently not used, no equatorial symmetries assumed, 
!            only whole sphere computations.    
!
!     mtrunc the highest longitudinal wave number retained in the
!            projection. It must be less than or equal to
!            the minimum of nlat-1 and nlon/2. The first wave
!            number is zero. For example, if wave numbers 0 and
!            1 are desired then mtrunc = 1.
!
!     lwshp  the dimension of the array wshp as it appears in the
!            program that calls shpgi. It must be at least
!            2*(nlat+1)**2+nlon+log2(nlon)
!
!     liwshp the dimension of the array iwshp as it appears in the
!            program that calls shpgi. It must be at least
!            4*(nlat+1).
!
!     lwork  the dimension of the array work as it appears in the
!            program that calls shpgi. It must be at least
!            1.25*(nlat+1)**2+7*nlat+8.
!
!     **************************************************************
!
!     output parameters
!
!     wshp   a single precision array that must be saved for
!            repeated use by subroutine shpg.        
!
!     iwshp  an integer array that must be saved for repeated
!            use by subroutine shpg.        
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
module module_shpg

    use spherepack_precision, only: &
        wp, & ! working precision
        ip, & ! integer precision
        PI

    use type_HFFTpack, only: &
    HFFTpack

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: shpg
    public :: shpgi

contains

subroutine shpgi(nlat, nlon, isym, mtrunc, wshp, lwshp, iwshp, &
 liwshp, work, lwork, ierror)


    type(HFFTpack) :: hfft
integer :: ierror
integer :: isym
integer :: iw1
integer :: iw2
integer :: iw3
integer :: iw4
integer :: iwshp
integer :: jw1
integer :: jw2
integer :: jw3
integer :: jw4
integer :: ktot
integer :: kw1
integer :: kw10
integer :: kw11
integer :: kw2
integer :: kw3
integer :: kw4
integer :: kw5
integer :: kw6
integer :: kw7
integer :: kw8
integer :: kw9
integer :: liwshp
integer :: log2n
integer :: lw1
integer :: lwork
integer :: lwshp
integer :: mlwk
integer :: mmax
integer :: mtrunc
integer :: nlat
integer :: nloc1
integer :: nloc2
integer :: nlon
integer :: nte
real :: wshp
real work(*)
dimension wshp(*), iwshp(*)
!
ierror = 1
if (nlat<1) return
ierror = 2
if (nlon<1) return
!      ierror = 3
!      if (isym.lt.0_wp .or. isym.gt.2) return
ierror = 4
mmax = min(nlat-1, nlon/2)
if (mtrunc<0 .or. mtrunc>mmax) return
ierror = 5
lw1 = 2*(nlat+1)**2
log2n = log(real(nlon))/log(2.0_wp)
if (lwshp<lw1+nlon+log2n) return
ierror = 6
if (liwshp<4*(nlat+1)) return
ierror = 7
mlwk = 1.25*(nlat+1)**2+7*nlat+8
if (lwork <mlwk) return
ierror = 0
!
call hfft%initialize(nlon, wshp(lw1+1))
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
kw4 = kw3+2*nte
kw5 = kw4+2*nte
kw6 = kw5+nte
kw7 = kw6+nte
kw8 = kw7+4*nte
kw9 = kw8+2*nte
kw10 = kw9+nloc1
kw11 = kw10+nloc1
ktot = kw11+nte*nte
!
call shpgi1(nlat, nlon, isym, mtrunc, nte, ierror, wshp(iw1), wshp(iw2), &
  wshp(iw3), wshp(iw4), iwshp(jw1), iwshp(jw2), iwshp(jw3), &
  iwshp(jw4), work(kw1), work(kw2), work(kw3), work(kw4), work(kw5), &
  work(kw6), work(kw7), work(kw8), work(kw9), work(kw10), work(kw11))
return
end subroutine shpgi
subroutine shpgi1(nlat, nlon, isym, mtrunc, idp, ierror, &
  pe, po, ze, zo, ipse, jzse, ipso, jzso, &
  cp, wx, thet, gwts, xx, z, a, b, ped, pod, u)

real :: dfn
real :: dmax
integer :: i
integer :: idp
integer :: ierr
integer :: ierror
integer :: ip
integer :: ipse
integer :: ipso
integer :: isym
integer :: it
integer :: j
integer :: js
integer :: jzse
integer :: jzso
integer :: k
integer :: lock
integer :: lwork
integer :: m
integer :: modn
integer :: mp1
integer :: ms2
integer :: mtrunc
integer :: mxtr
integer :: n
integer :: nec
integer :: nem
integer :: nlat
integer :: nlon
integer :: nmx
integer :: noc
integer :: nom
integer :: ns2
integer :: nshe
integer :: nsho
integer :: nte
integer :: nto
real :: pe
real :: po
real :: sum1
real :: toe
real :: tusl
real :: ze
real :: zo
real :: zort
!
real summation, eps, a1, b1, c1, work
parameter (eps=epsilon(1.0_wp))
real cp(idp), wx(idp), &
  thet(nlat), gwts(nlat), xx(idp), z(idp), a(4*idp), &
  b(2*idp), ped(idp, idp, 2), pod(idp, idp, 2), u(idp, idp)
!
dimension pe(idp, idp, 2), po(idp, idp, 2), ze(idp, idp, 2), &
  zo(idp, idp, 2), &
  ipse(idp, 2), jzse(idp, 2), ipso(idp, 2), jzso(idp, 2), &
  nshe(2), nsho(2)
dimension zort(64, 64, 2)
!
ns2 = nlat/2
modn = nlat-ns2-ns2
nte = (nlat+1)/2
nto = nlat-nte
tusl = 0.
toe = 0.
!
!     compute gauss grid distribution
!
lwork = nlat+1
call compute_gaussian_latitudes_and_weightsp(nlat, thet, gwts, work, lwork, ierr)
if (ierr /= 0) write(*, 160) ierr
160 format(' error in compute_gaussian_latitudes_and_weights =', i5)
do i=1, nto
gwts(i) = gwts(i)+gwts(i)
end do
!
!     compute n**2 basis (even functions)
!
do n=1, nlat+nlat-2
dfn = n
a(n) = sqrt(dfn*(dfn+1.0_wp))
end do
do n=1, nlat-1
dfn = n
b(n) = sqrt((dfn+dfn+3.0_wp)/(dfn+dfn-1.0_wp))
end do
!
mxtr = min(nlat-1, nlon/2, mtrunc)
ip = 2
do 200 mp1=1, mxtr+1
m = mp1-1
ip = 3-ip
ms2 = mp1/2
nem = (nlat-m+1)/2
nec = nte-nem
!
!     compute associated legendre functions
!
if (m<=1) then
do 205 j=1, nem
n = 2*j+m-2
call dlfkg(m, n, cp)
do i=1, nte
call dlftg (m, n, thet(i), cp, ped(i, j+nec, ip))
end do
205 continue
!
else
!
do 207 j=1, nem
n = 2*j+m-2
if (m>1.and.n>mxtr) then
do i=1, nte
u(i, j+nec) = ped(i, j+nec, ip)
end do
goto 207
end if
a1 = b(n-1)*a(n+m-3)/a(n+m-1)
b1 = a(n-m+1)/a(n+m-1)
if (n-m<=1) then
do i=1, nte
u(i, j+nec) = a1*ped(i, j+nec-1, ip) &
                   - b1*ped(i, j+nec, ip)    
end do
else
c1 = b(n-1)*a(n-m-1)/a(n+m-1)
do i=1, nte
u(i, j+nec) = a1*ped(i, j+nec-1, ip) &
   - b1*ped(i, j+nec, ip) + c1*u(i, j+nec-1)    
end do
end if
207 continue
do j=1, nem
do i=1, nte
ped(i, j+nec, ip) = u(i, j+nec)
end do
end do
end if
if (nec<=0) goto 200
!
!     generate orthogonal vector with 
!     random numbers using Fortran90
!     intrinsics RANDOM_{SEED, NUMBER}
!
!     comment out old code
!
!     do i=1, nte
!     xx(i) = rand()
!     end do
!
! replacement code
!
call RANDOM_SEED()
call random_number(xx(1:nte))
!
it = 0
201 do i=1, nte
z(i) = 0.0_wp
wx(i) = gwts(i)*xx(i)
end do
do 220 j=1, nte
if (j==nec) goto 220
call gs(nte, wx, ped(1, j, ip), z)
220 continue
!  
do i=1, nte
xx(i) = xx(i)-z(i)
end do
call normal(nte, xx, idp, gwts)
it = it+1
if (it<=2) goto 201
do i=1, nte
ped(i, nec, ip) = xx(i)
end do
200 continue
!
!     reorder if mtrunc is less than nlat-1 
!         case of even functions
!
nmx = nlat-mxtr
if (modn==1) then
nshe(1) = nmx/2
nshe(2) = (nmx-1)/2
else
nshe(1) = (nmx-1)/2
nshe(2) = nmx/2
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
summation = ped(j, i, ip)*gwts(j)
ze(j, i, ip) =  summation
pe(i, j, ip) = ped(i, j, ip)
if (abs(summation)>eps .and. lock==0) then
lock = 1
jzse(i, ip) = j
end if
end do
end do
250 continue
!
!     check orthogonality of pe(i, j, mp1)  mp1=1, 2
!
do ip=1, 2
dmax = 0.
do i=1, nte
do j=1, nte
sum1 = 0.
do k=1, nte
sum1 = sum1+ze(k, i, ip)*pe(k, j, ip)
end do
zo(i, j, ip) = sum1
if (i/=j) then
dmax = max(dmax, abs(sum1))
else
dmax = max(dmax, abs(sum1-1.0_wp))
end if
end do
end do
end do
!
!     compute n**2 basis (odd functions)
!
ip = 2
do 300 mp1=1, mxtr+1
ip = 3-ip
m = mp1-1
ms2 = mp1/2
nem = (nlat-m+1)/2
nom = nlat-m-nem 
noc = nto-nom
!
!     compute associated legendre functions
!
if (m<=1) then
do 305 j=1, nom
n = 2*j+m-1
call dlfkg(m, n, cp)
do i=1, nte
call dlftg (m, n, thet(i), cp, pod(i, j+noc, ip))
end do
if (modn>0) pod(nte, j+noc, ip) = 0.0_wp
305 continue
!
else
!
do 307 j=1, nom
n = 2*j+m-1
if (m>1.and.n>mxtr) then
do i=1, nte
u(i, j+noc) = pod(i, j+noc, ip)
end do
goto 304
end if
a1 = b(n-1)*a(n+m-3)/a(n+m-1)
b1 = a(n-m+1)/a(n+m-1)
if (n-m<=1) then
do i=1, nte
u(i, j+noc) = a1*pod(i, j+noc-1, ip) &
                   - b1*pod(i, j+noc, ip)    
end do
else
c1 = b(n-1)*a(n-m-1)/a(n+m-1)
do i=1, nte
u(i, j+noc) = a1*pod(i, j+noc-1, ip) &
   - b1*pod(i, j+noc, ip) + c1*u(i, j+noc-1)    
end do
end if
304 if (modn==1) u(nte, j+noc) = 0.0_wp
307 continue
do j=1, nom
do i=1, nte
pod(i, j+noc, ip) = u(i, j+noc)
end do
end do
end if
!
if (noc<=0) goto 300
!
!     old code with nonstandard (s)rand
!     commented out
!
!     do i=1, nte
!     xx(i) = rand()
!     end do
!
!     replacement code with standard Fortran90 
!     intrinsic
!
call random_number(xx(1:nte))
if (modn==1) xx(nte) = 0.0_wp
it = 0
306 do i=1, nte
z(i) = 0.
wx(i) = gwts(i)*xx(i)
end do
do 330 j=1, nto
if (j==noc) goto 330
call gs(nte, wx, pod(1, j, ip), z(1))
330 continue
!  
do i=1, nte
xx(i) = xx(i)-z(i)
end do
call normal(nte, xx, idp, gwts)
it = it+1
if (it<=2) goto 306
do i=1, nte
pod(i, noc, ip) = xx(i)
end do
if (modn==1) pod(nte, noc, ip) = 0.0_wp
300 continue
!
nmx = nlat-mxtr
if (modn==1) then
nsho(1) = (nmx-1)/2
nsho(2) = nmx/2
else
nsho(1) = nmx/2
nsho(2) = (nmx-1)/2
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
summation = pod(j, i, ip)*gwts(j)
zo(j, i, ip) = summation
po(i, j, ip) = pod(i, j, ip)
if (abs(summation)>eps .and. lock==0) then
lock = 1
jzso(i, ip) = j
end if
end do
end do
end do
!
!     check orthogonality of po(i, j, mp1)  mp1=1, 2
!
do ip=1, 2
dmax = 0.
do i=1, nto
do j=1, nto
sum1 = 0.
do k=1, nto
sum1 = sum1+zo(k, i, ip)*po(k, j, ip)
end do
zort(i, j, ip) = sum1
if (i/=j) then
dmax = max(dmax, abs(sum1))
else
dmax = max(dmax, abs(sum1-1.0_wp))
end if
end do
end do
end do
return
end subroutine shpgi1
!
!
! ... file shpg.f
!
! ... files which must be loaded with shpg.f
!
!     type_HFFTpack.f
!
!     shpg computes the harmonic projection, which is
!     equivalent to a harmonic analysis (forward) followed
!     by a harmonic synthesis (backward transform).
!     shpg uses the n**2 projection or complement when appropriate
!     as well as  odd/even factorization and zero truncation on an
!     on a Gaussian distributed grid as defined in the JCP paper
!     "Generalized discrete spherical harmonic transforms" 
!     by Paul N. Swarztrauber and William F. Spotz
!     J. Comp. Phys., 159(2000) pp. 213-230.
!
!     subroutine shpg(nlat, nlon, isym, mtrunc, x, y, idxy, 
!    1        wshp, lwshp, iwshp, liwshp, work, lwork, ierror)
!
!     shpg projects the array x onto the set of functions represented
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
!            nlon must be at least 4. 
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
!            appear in the program that calls shpg. It must be
!            at least nlat. 
!
!     wshp   a single precision array that must be saved for
!            repeated use by subroutine shpg.        
!
!     lwshp  the dimension of the array wshp as it appears in the
!            program that calls shpgi. It must be at least
!            2*(nlat+1)**2+nlon+log2(nlon)
!
!     iwshp  an integer array that must be saved for repeated
!            use by subroutine shpg.        
!
!
!     liwshp the dimension of the array iwshp as it appears in the
!            program that calls shpgi. It must be at least
!            4*(nlat+1).
!
!     work   a single precision work array that does 
!            not have to be saved.
!
!     lwork  the dimension of the array work as it appears in the
!            program that calls shpg. It must be at least
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
subroutine shpg(nlat, nlon, isym, mtrunc, x, y, idxy, &
        wshp, lwshp, iwshp, liwshp, work, lwork, ierror)


    type(HFFTpack) :: hfft
integer :: i
integer :: idxy
integer :: ierror
integer :: isym
integer :: iw1
integer :: iw2
integer :: iw3
integer :: iw4
integer :: iwshp
integer :: j
integer :: jw1
integer :: jw2
integer :: jw3
integer :: jw4
integer :: liwshp
integer :: log2n
integer :: lw1
integer :: lwork
integer :: lwshp
integer :: mmax
integer :: mtrunc
integer :: mwrk
integer :: nlat
integer :: nloc1
integer :: nloc2
integer :: nlon
integer :: nte
real :: sn
real :: work
real :: wshp
real :: x
real :: y
!
dimension wshp(*), iwshp(*), work(*), x(idxy, nlon), y(idxy, nlon)
!
ierror = 1
if (nlat<1) return
ierror = 2
if (nlon<1) return
!      ierror = 3
!      if (isym.lt.0_wp .or. isym.gt.2) return
ierror = 4
mmax = min(nlat-1, nlon/2)
if (mtrunc<0 .or. mtrunc>mmax) return
ierror = 5
log2n = log(real(nlon))/log(2.0_wp)
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
call hfft%forward(nlat, nlon, y, idxy, wshp(lw1+1), work)
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
call shpg1(nlat, nlon, isym, mtrunc, y, y, idxy, ierror, &
 nte, wshp(iw1), wshp(iw2), wshp(iw3), wshp(iw4), iwshp(jw1), &
 iwshp(jw2), iwshp(jw3), iwshp(jw4), work(jw1), &
 work(jw2), work(jw3), work(jw4))
!
call hfft%backward(nlat, nlon, y, idxy, wshp(lw1+1), work)
!
sn = 1.0_wp/nlon
do j=1, nlon
 do i=1, nlat
  y(i, j) = sn*y(i, j)
 end do
end do
return
end subroutine shpg
subroutine shpg1(nlat, nlon, isym, mtrunc, sx, sy, idxy, ierror, &
 idp, pe, po, ze, zo, ipse, jzse, ipso, jzso, xe, xo, ye, yo)

integer :: i
integer :: idp
integer :: idxy
integer :: ierror
integer :: ip
integer :: ipse
integer :: ipso
integer :: isym
integer :: j
integer :: js
integer :: jzse
integer :: jzso
integer :: lag
integer :: m
integer :: modn
integer :: mp1
integer :: mpm
integer :: ms2
integer :: mtrunc
integer :: mxtr
integer :: nec
integer :: nem
integer :: nlat
integer :: nlon
integer :: nmx
integer :: noc
integer :: nom
integer :: ns2
integer :: nshe
integer :: nsho
integer :: nte
integer :: nto
real :: pe
real :: po
real :: sx
real :: sy
real :: xe
real :: xo
real :: ye
real :: yo
real :: ze
real :: zo
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
mxtr = min(nlat-1, nlon/2, mtrunc)
nmx = nlat-mxtr
if (modn==1) then
nshe(1) = nmx/2
nshe(2) = (nmx-1)/2
nsho(1) = (nmx-1)/2
nsho(2) = nmx/2
else
nshe(1) = (nmx-1)/2
nshe(2) = nmx/2
nsho(1) = nmx/2
nsho(2) = (nmx-1)/2
end if
!
ip = 2 
do 100 mp1=1, mxtr+1
ip = 3-ip
if (mxtr==nlat-1.and.mp1==1) then
do i=1, nlat
sy(i, mp1) = sx(i, mp1)
end do
!      if (mp1.eq.2) then
!      sy(1, 2) = 0.
!      sy(nlat, 2) = 0.
!      end if
!      if (nlon.ge.3) then
!      sy(1, 3) = 0.
!      sy(nlat, 3) = 0.
!      do i=2, nlat-1
!      sy(i, 3) = sx(i, 3)
!      end do
!      end if
goto 100
end if
m = mp1-1
mpm = max(1, m+m)
ms2 = mp1/2
!      mrank = min(nlat-m, nlat-ms2-ms2)   
!      nrank = nlat-mrank
!      nem = (mrank+1)/2-nshe(ip)
!      nom = mrank-(mrank+1)/2-nsho(ip)
nem = (nlat-m+1)/2-nshe(ip)
nom = (nlat-m)/2-nsho(ip)
nec = nte-nem
noc = nto-nom
do i=1, nte
xe(i, 1) = .5*(sx(i, mpm)+sx(nlat+1-i, mpm))
xo(i, 1) = .5*(sx(i, mpm)-sx(nlat+1-i, mpm))
end do
!      if (modn.eq.1) then
!      xe(nte, 1) = sx(nte, mpm)
!      xo(nte, 1) = 0.
!      end if
if (mpm<nlon) then
do i=1, nte
xe(i, 2) = .5*(sx(i, mpm+1)+sx(nlat+1-i, mpm+1))
xo(i, 2) = .5*(sx(i, mpm+1)-sx(nlat+1-i, mpm+1))
end do
!      if (modn.eq.1) then
!      xe(nte, 2) = sx(nte, mpm+1)
!      xo(nte, 2) = 0.
!      end if
end if
lag = 0
if (m==0.or.mpm==nlon) lag = 1
if (3*nec<2*nem.or.nem==0) then
call tmxmx(lag, nte, nec, idp, pe(1, 1, ip), nte, idp, &
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
call tmxmx(lag, nte, nem, idp, pe(1, nec+1, ip), nte, idp, &
ze(1, nec+1, ip), xe, ye, ipse(nec+1, ip), jzse(nec+1, ip))
end if
if (3*noc<2*nom.or.nom==0) then
call tmxmx(lag, nto, noc, idp, po(1, 1, ip), nto, idp, &
          zo(1, 1, ip), xo, yo, ipso(1, ip), jzso(1, ip))
do i=1, nto
yo(i, 1) = xo(i, 1)-yo(i, 1)
end do
if (mpm<nlon.and.m/=0) then
do i=1, nto
yo(i, 2) = xo(i, 2)-yo(i, 2)
end do
end if
else
call tmxmx(lag, nto, nom, idp, po(1, noc+1, ip), nto, idp, &
zo(1, noc+1, ip), xo, yo, ipso(noc+1, ip), jzso(noc+1, ip))  
end if
do i=1, nto
sy(i, mpm) = ye(i, 1)+yo(i, 1)
sy(nlat+1-i, mpm) = ye(i, 1)-yo(i, 1)
end do
if (nte>nto) sy(nte, mpm) = ye(nte, 1)
if (mpm<nlon.and.m/=0) then
do i=1, nto
sy(i, mpm+1) = ye(i, 2)+yo(i, 2)
sy(nlat+1-i, mpm+1) = ye(i, 2)-yo(i, 2)
end do 
if (nte>nto) sy(nte, mpm+1) = ye(nte, 2)
end if
100 continue
!
js = mxtr+mxtr+2
do j=js, nlon
do i=1, nlat
sy(i, j) = 0.
end do
end do
return
end subroutine shpg1
subroutine mxm(lr, lc, ld, a, mc, md, b, nd, c)

integer :: i
integer :: j
integer :: k
integer :: lc
integer :: ld
integer :: lr
integer :: mc
integer :: md
integer :: nd
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

real :: a
real :: b
real :: c
integer :: i
integer :: j
integer :: k
integer :: lc
integer :: ld
integer :: lr
integer :: mc
integer :: md
integer :: nd
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

real :: a
real :: b
integer :: i
integer :: j
integer :: k
integer :: lc
integer :: ld
integer :: lr
integer :: mc
integer :: md
real :: sum1
real :: sum2
real :: x
real :: y
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

integer :: i
integer :: j
integer :: k
integer :: lc
integer :: ld
integer :: lr
integer :: mc
integer :: md
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
subroutine tmxmx(lag, lr, lc, ld, a, mc, md, b, x, y, is, js)

real :: a
real :: b
integer :: i
integer :: is
integer :: j
integer :: js
integer :: k
integer :: kmx
integer :: lag
integer :: lc
integer :: ld
integer :: lr
integer :: mc
integer :: md
real :: sum1
real :: sum2
real :: x
real :: y
dimension a(ld, *), b(md, *), x(ld, 2), y(ld, 2), &
              is(*), js(*)
!
kmx = min(lr+1, ld)
if (lag==1) then
do k=1, kmx
y(k, 1) = 0.
end do
!      if (lc.eq.0_wp) then
!      do k=1, lr
!      y(k, 1) = x(k, 1)
!      end do
!      return
!      end if
if (lc<=0) return
do i=1, lc
sum1 = 0.
do j=js(i), mc
sum1 = sum1 + b(j, i)*x(j, 1)
end do
do k=is(i), lr
y(k, 1) = y(k, 1)+sum1*a(k, i)
end do
end do
return
end if
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

integer :: i
integer :: idp
integer :: ijs
integer :: irc
integer :: j
integer :: n
integer :: nrc
real a, eps
parameter (eps=epsilon(1.0_wp))
dimension a(idp, *), ijs(n)
!
!     irc = 0 for columns , or irc = 1 for rows
!
if (irc/=0) goto 30
do 20 j=1, nrc
do i=1, n
ijs(j) = i
if (abs(a(i, j)) > eps) goto 20
end do
20 continue
return
30 do 50 i=1, nrc
do j=1, n
ijs(i) = j
if (abs(a(i, j)) > eps) goto 50
end do
50 continue
return
end subroutine trunc
subroutine gs(n, x, y, z)

integer :: i
integer :: n
dimension x(n), y(n), z(n)
real x, y, z, summation
!
!     accumulate innerproducts of x with respect to y.
!
summation = 0.
do i=1, n
summation = summation+x(i)*y(i)
end do
do i=1, n
z(i) = z(i)+summation*y(i)
end do
return
end subroutine gs
subroutine normal(n, x, id, q)

integer :: i
integer :: id
integer :: n
dimension x(n), q(n)
real x, q, sqs
!
!     normalize x
!
sqs = 0.
do i=1, n
!      sum = 0.
!      do j=1, n
!      sum = sum+q(i, j)*x(j)
!      end do
!      sqs = sqs+sum*x(i)
sqs = sqs+q(i)*x(i)*x(i)
end do
!
if (sqs /= 0) goto 4
write(*, 3)
3 format(' norm of z is zero in subroutine normal')
return
4 sqs = sqrt(sqs)
do i=1, n
x(i) = x(i)/sqs
end do
return
end subroutine normal
subroutine coe(moe, n, x, dmax)

integer :: i
integer :: moe
integer :: n
integer :: nh
real x(n), dmax
nh = (n+1)/2
dmax = 0.
if (moe/=0) goto 1
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
!     subroutine dlfkg(m, n, cp)
!
!     subroutine dlfkg computes the coefficients in the trigonometric
!     expansion of the normalized associated legendre functions:
!
!     pbar(m, n, theta) = sqrt((2*n+1)*factorial(n-m)/(2*factorial(n+m)))
!                        *sin(theta)**m/(2**n*factorial(n)) times the
!                        (n+m)th derivative of (x**2-1)**n with respect
!                        to x=cos(theta)
!
!     where theta is colatitude.
!
!     subroutine dlfkg computes the coefficients cp(k) in the
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
subroutine dlfkg (m, n, cp)

integer :: i
integer :: l
integer :: m
integer :: ma
integer :: n
integer :: nex
integer :: nmms2
!
real cp, fnum, fden, fnmh, a1, b1, c1, cp2, fnnp1, fnmsq, fk, &
       t1, t2, pm1, sc10, sc20, sc40
dimension       cp(1)
parameter (sc10=1024.0_wp)
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
2 cp(1) = sqrt(2.0_wp)
return
3 if (ma /= 0) goto 4
cp(1) = sqrt(1.5)
return
4 cp(1) = sqrt(.75)
if (m == -1) cp(1) = -cp(1)
return
5 if (mod(n+ma, 2) /= 0) goto 10
nmms2 = (n-ma)/2
fnum = n+ma+1
fnmh = n-ma+1
pm1 = 1.0_wp
goto 15
10 nmms2 = (n-ma-1)/2
fnum = n+ma+2
fnmh = n-ma+2
pm1 = -1.0_wp
15 t1 = 1.0_wp/sc20
nex = 20
fden = 2.0_wp
if (nmms2 < 1) goto 20
do 18 i=1, nmms2
t1 = fnum*t1/fden
if (t1 > sc20) then
t1 = t1/sc40
nex = nex+40
end if
fnum = fnum+2.
fden = fden+2.
18 continue
20 t1 = t1/2.0_wp**(n-1-nex)
if (mod(ma/2, 2) /= 0) t1 = -t1
t2 = 1. 
if (ma == 0) goto 26
do 25 i=1, ma
t2 = fnmh*t2/(fnmh+pm1)
fnmh = fnmh+2.
25 continue
26 cp2 = t1*sqrt((n+.5)*t2)
fnnp1 = n*(n+1)
fnmsq = fnnp1-2.0_wp*ma*ma
l = (n+1)/2
if (mod(n, 2) == 0 .and. mod(ma, 2) == 0) l = l+1
cp(l) = cp2
if (m >= 0) goto 29
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
goto 30
end subroutine dlfkg
subroutine dlftg (m, n, theta, cp, pb)

integer :: k
integer :: kdo
integer :: m
integer :: mmod
integer :: n
integer :: nmod
dimension cp(1)
real cp, pb, theta, cdt, sdt, cth, sth, chh
cdt = cos(2.0_wp*theta)
sdt = sin(2.0_wp*theta)
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
end subroutine dlftg
!
subroutine compute_gaussian_latitudes_and_weightsp(nlat, theta, wts, w, lwork, ierror)

real :: eps
integer :: i
integer :: idx
integer :: ierror
integer :: it
integer :: lwork
integer :: mnlat
integer :: nhalf
integer :: nix
integer :: nlat
integer :: ns2
real :: sgnd
!
!                             April 2002
!                        
!     gauss points and weights are computed using the fourier-newton
!     described in "on computing the points and weights for 
!     gauss-legendre quadrature", paul n. swarztrauber, siam journal 
!     on scientific computing that has been accepted for publication.
!     This routine is faster and more accurate than older program 
!     with the same name.
!
!     subroutine compute_gaussian_latitudes_and_weightsp computes the nlat gaussian colatitudes and weights
!     in real. the colatitudes are in radians and lie in the
!     in the interval (0, pi).
!
!     input parameters
!
!     nlat    the number of gaussian colatitudes in the interval (0, pi)
!             (between the two poles).  nlat must be greater than zero.
!
!     w       unused variable that permits a simple exchange with the
!             old routine with the same name in spherepack.
!
!     lwork   unused variable that permits a simple exchange with the
!             old routine with the same name in spherepack.
!
!     output parameters
!
!     theta   a real array with length nlat
!             containing the gaussian colatitudes in
!             increasing radians on the interval (0, pi).
!
!     wts     a real array with lenght nlat
!             containing the gaussian weights.
!
!     ierror = 0 no errors
!            = 1 if nlat.le.0_wp
!
!  *****************************************************************
!
real theta(nlat), wts(nlat), &
 x, HALF_PI, dtheta, dthalf, cmax, zprev, zlast, zero, &
 zhold, pb, dpb, dcor, summation, w, cz
!
!     check work space length
!
ierror = 1
if (nlat<=0) return
ierror = 0
!
!     compute weights and points analytically when nlat=1, 2
!
if (nlat==1) then
theta(1) = PI/2
wts(1) = 2.0_wp
return
end if
if (nlat==2) then
x = sqrt(1.0_wp/3.0_wp)
theta(1) = acos(x)
theta(2) = acos(-x)
wts(1) = 1.0_wp
wts(2) = 1.0_wp
return
end if
eps = sqrt(epsilon(1.0_wp))
eps = eps*sqrt(eps)
HALF_PI = pi/2
mnlat = mod(nlat, 2)
ns2 = nlat/2
nhalf = (nlat+1)/2
idx = ns2+2
!
call cpdp1 (nlat, cz, theta(ns2+1), wts(ns2+1))
!
dtheta = HALF_PI/nhalf
dthalf = dtheta/2.0_wp
cmax = .2*dtheta
!
!     estimate first point next to theta = pi/2
!
if (mnlat/=0) then
zero = HALF_PI-dtheta
zprev = HALF_PI
nix = nhalf-1
else
zero = HALF_PI-dthalf
nix = nhalf
end if
9 it = 0
10 it = it+1
zlast = zero
!
!     newton iterations
!
call tpdp1 (nlat, zero, cz, theta(ns2+1), wts(ns2+1), pb, dpb)
dcor = pb/dpb
sgnd = 1.0_wp
if (dcor /= 0.0_wp) sgnd = dcor/abs(dcor)
dcor = sgnd*min(abs(dcor), cmax)
zero = zero-dcor
if (abs(zero-zlast)>eps*abs(zero)) goto 10
theta(nix) = zero
zhold = zero
!      wts(nix) = (nlat+nlat+1)/(dpb*dpb)
!    
!     yakimiw's formula permits using old pb and dpb
!
wts(nix) = (nlat+nlat+1)/(dpb+pb*cos(zlast)/sin(zlast))**2
nix = nix-1
if (nix==0) goto 30
if (nix==nhalf-1)  zero = 3.0_wp*zero-pi
if (nix<nhalf-1)  zero = zero+zero-zprev
zprev = zhold
goto 9
!
!     extend points and weights via symmetries
!
30 if (mnlat/=0) then
theta(nhalf) = HALF_PI
call tpdp1 (nlat, HALF_PI, cz, theta(ns2+1), wts(ns2+1), pb, dpb)
wts(nhalf) = (nlat+nlat+1)/(dpb*dpb)
end if
do i=1, ns2
wts(nlat-i+1) = wts(i)
theta(nlat-i+1) = pi-theta(i)
end do
summation = 0.0_wp
do i=1, nlat
summation = summation+wts(i)
end do
do i=1, nlat
wts(i) = 2.0_wp*wts(i)/summation
end do
return
end subroutine compute_gaussian_latitudes_and_weightsp
subroutine cpdp1(n, cz, cp, dcp)

integer :: j
integer :: n
integer :: ncp
!
!     computes the fourier coefficients of the legendre
!     polynomial p_n^0 and its derivative. 
!     n is the degree and n/2 or (n+1)/2
!     coefficients are returned in cp depending on whether
!     n is even or odd. The same number of coefficients
!     are returned in dcp. For n even the constant 
!     coefficient is returned in cz. 
!
real cp(n/2+1), dcp(n/2+1), &
 t1, t2, t3, t4, cz
ncp = (n+1)/2
t1 = -1.0_wp
t2 = n+1.0_wp
t3 = 0.0_wp
t4 = n+n+1.0_wp
if (mod(n, 2)==0) then
cp(ncp) = 1.0_wp
do j = ncp, 2, -1
t1 = t1+2.0_wp
t2 = t2-1.0_wp
t3 = t3+1.0_wp
t4 = t4-2.0_wp
cp(j-1) = (t1*t2)/(t3*t4)*cp(j)
end do
t1 = t1+2.0_wp
t2 = t2-1.0_wp
t3 = t3+1.0_wp
t4 = t4-2.0_wp
cz = (t1*t2)/(t3*t4)*cp(1)
do j=1, ncp
dcp(j) = (j+j)*cp(j)
end do
else
cp(ncp) = 1.0_wp
do j = ncp-1, 1, -1
t1 = t1+2.0_wp
t2 = t2-1.0_wp
t3 = t3+1.0_wp
t4 = t4-2.0_wp
cp(j) = (t1*t2)/(t3*t4)*cp(j+1)
end do
do j=1, ncp
dcp(j) = (j+j-1)*cp(j)
end do
end if
return
end subroutine cpdp1
subroutine tpdp1 (n, theta, cz, cp, dcp, pb, dpb)

integer :: k
integer :: kdo
integer :: n
!
!     computes pn(theta) and its derivative dpb(theta) with 
!     respect to theta
!
real cp(n/2+1), dcp(n/2+1), cz, &
  pb, dpb, fn, theta, cdt, sdt, cth, sth, chh
!
fn = n
cdt = cos(2.0_wp*theta)
sdt = sin(2.0_wp*theta)
if (mod(n, 2) ==0) then
!
!     n even
!
kdo = n/2
pb = .5*cz
dpb = 0.0_wp
if (n > 0) then
cth = cdt
sth = sdt
do 170 k=1, kdo
!      pb = pb+cp(k)*cos(2*k*theta)
pb = pb+cp(k)*cth
!      dpb = dpb-(k+k)*cp(k)*sin(2*k*theta)
dpb = dpb-dcp(k)*sth
chh = cdt*cth-sdt*sth
sth = sdt*cth+cdt*sth
cth = chh
170 continue
end if
else
!
!     n odd
!
kdo = (n+1)/2
pb = 0.0_wp
dpb = 0.0_wp
cth = cos(theta)
sth = sin(theta)
do 190 k=1, kdo
!      pb = pb+cp(k)*cos((2*k-1)*theta)
pb = pb+cp(k)*cth
!      dpb = dpb-(k+k-1)*cp(k)*sin((2*k-1)*theta)
dpb = dpb-dcp(k)*sth
chh = cdt*cth-sdt*sth
sth = sdt*cth+cdt*sth
cth = chh
190 continue
end if
return
end subroutine tpdp1

end module module_shpg
