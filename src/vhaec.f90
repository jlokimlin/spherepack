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
!
!
! ... file vhaec.f
!
!     this file contains code and documentation for subroutines
!     vhaec and vhaeci
!
! ... files which must be loaded with vhaec.f
!
!     type_SpherepackAux.f, type_HFFTpack.f
!
!                                                                              
!     subroutine vhaec(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, 
!    +                 mdab, ndab, wvhaec, lvhaec, work, lwork, ierror)
!
!     subroutine vhaec performs the vector spherical harmonic analysis
!     on the vector field (v, w) and stores the result in the arrays
!     br, bi, cr, and ci. v(i, j) and w(i, j) are the colatitudinal 
!     (measured from the north pole) and east longitudinal components
!     respectively, located at colatitude theta(i) = (i-1)*pi/(nlat-1)
!     and longitude phi(j) = (j-1)*2*pi/nlon. the spectral
!     representation of (v, w) is given at output parameters v, w in 
!     subroutine vhsec.  
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
!            than zero. the axisymmetric case corresponds to nlon=1.
!            the efficiency of the computation is improved when nlon
!            is a product of small prime numbers.
!
!     ityp   = 0  no symmetries exist about the equator. the analysis
!                 is performed on the entire sphere.  i.e. on the
!                 arrays v(i, j), w(i, j) for i=1, ..., nlat and 
!                 j=1, ..., nlon.   
!
!            = 1  no symmetries exist about the equator. the analysis
!                 is performed on the entire sphere.  i.e. on the
!                 arrays v(i, j), w(i, j) for i=1, ..., nlat and 
!                 j=1, ..., nlon. the curl of (v, w) is zero. that is, 
!                 (d/dtheta (sin(theta) w) - dv/dphi)/sin(theta) = 0. 
!                 the coefficients cr and ci are zero.
!
!            = 2  no symmetries exist about the equator. the analysis
!                 is performed on the entire sphere.  i.e. on the
!                 arrays v(i, j), w(i, j) for i=1, ..., nlat and 
!                 j=1, ..., nlon. the divergence of (v, w) is zero. i.e., 
!                 (d/dtheta (sin(theta) v) + dw/dphi)/sin(theta) = 0. 
!                 the coefficients br and bi are zero.
!
!            = 3  v is symmetric and w is antisymmetric about the 
!                 equator. the analysis is performed on the northern
!                 hemisphere only.  i.e., if nlat is odd the analysis
!                 is performed on the arrays v(i, j), w(i, j) for 
!                 i=1, ..., (nlat+1)/2 and j=1, ..., nlon. if nlat is
!                 even the analysis is performed on the the arrays
!                 v(i, j), w(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
!
!            = 4  v is symmetric and w is antisymmetric about the 
!                 equator. the analysis is performed on the northern
!                 hemisphere only.  i.e., if nlat is odd the analysis
!                 is performed on the arrays v(i, j), w(i, j) for 
!                 i=1, ..., (nlat+1)/2 and j=1, ..., nlon. if nlat is
!                 even the analysis is performed on the the arrays
!                 v(i, j), w(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
!                 the curl of (v, w) is zero. that is, 
!                 (d/dtheta (sin(theta) w) - dv/dphi)/sin(theta) = 0. 
!                 the coefficients cr and ci are zero.
!
!            = 5  v is symmetric and w is antisymmetric about the 
!                 equator. the analysis is performed on the northern
!                 hemisphere only.  i.e., if nlat is odd the analysis
!                 is performed on the arrays v(i, j), w(i, j) for 
!                 i=1, ..., (nlat+1)/2 and j=1, ..., nlon. if nlat is
!                 even the analysis is performed on the the arrays
!                 v(i, j), w(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
!                 the divergence of (v, w) is zero. i.e., 
!                 (d/dtheta (sin(theta) v) + dw/dphi)/sin(theta) = 0. 
!                 the coefficients br and bi are zero.
!
!            = 6  v is antisymmetric and w is symmetric about the 
!                 equator. the analysis is performed on the northern
!                 hemisphere only.  i.e., if nlat is odd the analysis
!                 is performed on the arrays v(i, j), w(i, j) for 
!                 i=1, ..., (nlat+1)/2 and j=1, ..., nlon. if nlat is
!                 even the analysis is performed on the the arrays
!                 v(i, j), w(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
!
!            = 7  v is antisymmetric and w is symmetric about the 
!                 equator. the analysis is performed on the northern
!                 hemisphere only.  i.e., if nlat is odd the analysis
!                 is performed on the arrays v(i, j), w(i, j) for 
!                 i=1, ..., (nlat+1)/2 and j=1, ..., nlon. if nlat is
!                 even the analysis is performed on the the arrays
!                 v(i, j), w(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
!                 the curl of (v, w) is zero. that is, 
!                 (d/dtheta (sin(theta) w) - dv/dphi)/sin(theta) = 0. 
!                 the coefficients cr and ci are zero.
!
!            = 8  v is antisymmetric and w is symmetric about the 
!                 equator. the analysis is performed on the northern
!                 hemisphere only.  i.e., if nlat is odd the analysis
!                 is performed on the arrays v(i, j), w(i, j) for 
!                 i=1, ..., (nlat+1)/2 and j=1, ..., nlon. if nlat is
!                 even the analysis is performed on the the arrays
!                 v(i, j), w(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
!                 the divergence of (v, w) is zero. i.e., 
!                 (d/dtheta (sin(theta) v) + dw/dphi)/sin(theta) = 0. 
!                 the coefficients br and bi are zero.
!
!
!     nt     the number of analyses.  in the program that calls vhaec, 
!            the arrays v, w, br, bi, cr, and ci can be three dimensional
!            in which case multiple analyses will be performed.
!            the third index is the analysis index which assumes the 
!            values k=1, ..., nt.  for a single analysis set nt=1. the
!            discription of the remaining parameters is simplified
!            by assuming that nt=1 or that all the arrays are two
!            dimensional.
!
!     v, w    two or three dimensional arrays (see input parameter nt)
!            that contain the vector function to be analyzed.
!            v is the colatitudnal component and w is the east 
!            longitudinal component. v(i, j), w(i, j) contain the
!            components at colatitude theta(i) = (i-1)*pi/(nlat-1)
!            and longitude phi(j) = (j-1)*2*pi/nlon. the index ranges
!            are defined above at the input parameter ityp.
!
!     idvw   the first dimension of the arrays v, w as it appears in
!            the program that calls vhaec. if ityp .le. 2 then idvw
!            must be at least nlat.  if ityp .gt. 2 and nlat is
!            even then idvw must be at least nlat/2. if ityp .gt. 2
!            and nlat is odd then idvw must be at least (nlat+1)/2.
!
!     jdvw   the second dimension of the arrays v, w as it appears in
!            the program that calls vhaec. jdvw must be at least nlon.
!
!     mdab   the first dimension of the arrays br, bi, cr, and ci as it
!            appears in the program that calls vhaec. mdab must be at
!            least min(nlat, nlon/2) if nlon is even or at least
!            min(nlat, (nlon+1)/2) if nlon is odd.
!
!     ndab   the second dimension of the arrays br, bi, cr, and ci as it
!            appears in the program that calls vhaec. ndab must be at
!            least nlat.
!
!     wvhaec an array which must be initialized by subroutine vhaeci.
!            once initialized, wvhaec can be used repeatedly by vhaec
!            as long as nlon and nlat remain unchanged.  wvhaec must
!            not be altered between calls of vhaec.
!
!     lvhaec the dimension of the array wvhaec as it appears in the
!            program that calls vhaec. define
!
!               l1 = min(nlat, nlon/2) if nlon is even or
!               l1 = min(nlat, (nlon+1)/2) if nlon is odd
!
!            and
!
!               l2 = nlat/2        if nlat is even or
!               l2 = (nlat+1)/2    if nlat is odd
!
!            then lvhaec must be at least
!
!            4*nlat*l2+3*max(l1-2, 0)*(nlat+nlat-l1-1)+nlon+15
!
!
!     work   a work array that does not have to be saved.
!
!     lwork  the dimension of the array work as it appears in the
!            program that calls vhaec. define
!
!               l2 = nlat/2        if nlat is even or
!               l2 = (nlat+1)/2    if nlat is odd
!
!            if ityp .le. 2 then lwork must be at least
!
!                    nlat*(2*nt*nlon+max(6*l2, nlon))
!
!            if ityp .gt. 2 then lwork must be at least
!
!                    l2*(2*nt*nlon+max(6*nlat, nlon))
!
!     **************************************************************
!
!     output parameters
!
!     br, bi  two or three dimensional arrays (see input parameter nt)
!     cr, ci  that contain the vector spherical harmonic coefficients
!            in the spectral representation of v(i, j) and w(i, j) given 
!            in the discription of subroutine vhsec. br(mp1, np1), 
!            bi(mp1, np1), cr(mp1, np1), and ci(mp1, np1) are computed 
!            for mp1=1, ..., mmax and np1=mp1, ..., nlat except for np1=nlat
!            and odd mp1. mmax=min(nlat, nlon/2) if nlon is even or 
!            mmax=min(nlat, (nlon+1)/2) if nlon is odd. 
!      
!     ierror = 0  no errors
!            = 1  error in the specification of nlat
!            = 2  error in the specification of nlon
!            = 3  error in the specification of ityp
!            = 4  error in the specification of nt
!            = 5  error in the specification of idvw
!            = 6  error in the specification of jdvw
!            = 7  error in the specification of mdab
!            = 8  error in the specification of ndab
!            = 9  error in the specification of lvhaec
!            = 10 error in the specification of lwork
!
!
! *******************************************************************
!
!     subroutine vhaeci(nlat, nlon, wvhaec, lvhaec, dwork, ldwork, ierror)
!
!     subroutine vhaeci initializes the array wvhaec which can then be
!     used repeatedly by subroutine vhaec until nlat or nlon is changed.
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
!            than zero. the axisymmetric case corresponds to nlon=1.
!            the efficiency of the computation is improved when nlon
!            is a product of small prime numbers.
!
!     lvhaec the dimension of the array wvhaec as it appears in the
!            program that calls vhaec. define
!
!               l1 = min(nlat, nlon/2) if nlon is even or
!               l1 = min(nlat, (nlon+1)/2) if nlon is odd
!
!            and
!
!               l2 = nlat/2        if nlat is even or
!               l2 = (nlat+1)/2    if nlat is odd
!
!            then lvhaec must be at least
!
!            4*nlat*l2+3*max(l1-2, 0)*(nlat+nlat-l1-1)+nlon+15
!
!
!     dwork  a real work array that does not have to be saved.
!
!     ldwork the dimension of the array dwork as it appears in the
!            program that calls vhaec. ldwork must be at least
!            2*(nlat+2)
!
!
!     **************************************************************
!
!     output parameters
!
!     wvhaec an array which is initialized for use by subroutine vhaec.
!            once initialized, wvhaec can be used repeatedly by vhaec
!            as long as nlat or nlon remain unchanged.  wvhaec must not
!            be altered between calls of vhaec.
!
!
!     ierror = 0  no errors
!            = 1  error in the specification of nlat
!            = 2  error in the specification of nlon
!            = 3  error in the specification of lvhaec
!            = 4  error in the specification of ldwork
!
!
module module_vhaec

    use spherepack_precision, only: &
        wp, & ! working precision
        ip ! integer precision

    use type_HFFTpack, only: &
        HFFTpack

    use type_SpherepackAux, only: &
        SpherepackAux

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    public :: vhaec
    public :: vhaeci

contains

subroutine vhaec(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, &
           mdab, ndab, wvhaec, lvhaec, work, lwork, ierror)
implicit none
real(wp) :: bi
real(wp) :: br
real(wp) :: ci
real(wp) :: cr
integer(ip) :: idv
integer(ip) :: idvw
integer(ip) :: ierror
integer(ip) :: imid
integer(ip) :: ist
integer(ip) :: ityp
integer(ip) :: iw1
integer(ip) :: iw2
integer(ip) :: iw3
integer(ip) :: iw4
integer(ip) :: iw5
integer(ip) :: jdvw
integer(ip) :: jw1
integer(ip) :: jw2
integer(ip) :: labc
integer(ip) :: lnl
integer(ip) :: lvhaec
integer(ip) :: lwork
integer(ip) :: lwzvin
integer(ip) :: lzz1
integer(ip) :: mdab
integer(ip) :: mmax
integer(ip) :: ndab
integer(ip) :: nlat
integer(ip) :: nlon
integer(ip) :: nt
real(wp) :: v
real(wp) :: w
real(wp) :: work
real(wp) :: wvhaec
dimension v(idvw, jdvw, nt), w(idvw, jdvw, nt), br(mdab, ndab, nt), &
          bi(mdab, ndab, nt), cr(mdab, ndab, nt), ci(mdab, ndab, nt), &
          work(lwork), wvhaec(lvhaec)

ierror = 1
if (nlat < 3) return
ierror = 2
if (nlon < 1) return
ierror = 3
if (ityp<0 .or. ityp>8) return
ierror = 4
if (nt < 0) return
ierror = 5
imid = (nlat+1)/2
if ((ityp<=2 .and. idvw<nlat) .or. &
   (ityp>2 .and. idvw<imid)) return
ierror = 6
if (jdvw < nlon) return
ierror = 7
mmax = min(nlat, (nlon+1)/2)
if (mdab < mmax) return
ierror = 8
if (ndab < nlat) return
ierror = 9
lzz1 = 2*nlat*imid
labc = 3*(max(mmax-2, 0)*(nlat+nlat-mmax-1))/2
if (lvhaec < 2*(lzz1+labc)+nlon+15) return
ierror = 10
if (ityp <= 2 .and. &
         lwork < nlat*(2*nt*nlon+max(6*imid, nlon))) return
if (ityp > 2 .and. &
         lwork < imid*(2*nt*nlon+max(6*nlat, nlon))) return
ierror = 0
idv = nlat
if (ityp > 2) idv = imid
lnl = nt*idv*nlon
ist = 0
if (ityp <= 2) ist = imid
iw1 = ist+1
iw2 = lnl+1
iw3 = iw2+ist
iw4 = iw2+lnl
iw5 = iw4+3*imid*nlat
lwzvin = lzz1+labc
jw1 = lwzvin+1
jw2 = jw1+lwzvin
call vhaec1(nlat, nlon, ityp, nt, imid, idvw, jdvw, v, w, mdab, ndab, &
     br, bi, cr, ci, idv, work, work(iw1), work(iw2), work(iw3), &
     work(iw4), work(iw5), wvhaec, wvhaec(jw1), wvhaec(jw2))

contains

subroutine vhaec1(nlat, nlon, ityp, nt, imid, idvw, jdvw, v, w, mdab, &
   ndab, br, bi, cr, ci, idv, ve, vo, we, wo, zv, zw, wzvin, wzwin, wrfft)

real(wp) :: bi
real(wp) :: br
real(wp) :: ci
real(wp) :: cr
real(wp) :: fsn
integer(ip) :: i
integer(ip) :: idv
integer(ip) :: idvw
integer(ip) :: imid
integer(ip) :: imm1
integer(ip) :: ityp

integer(ip) :: iv
integer(ip) :: iw
integer(ip) :: j
integer(ip) :: jdvw
integer(ip) :: k
integer(ip) :: m
integer(ip) :: mdab
integer(ip) :: mlat
integer(ip) :: mlon
integer(ip) :: mmax
integer(ip) :: mp1
integer(ip) :: mp2
integer(ip) :: ndab
integer(ip) :: ndo1
integer(ip) :: ndo2
integer(ip) :: nlat
integer(ip) :: nlon
integer(ip) :: nlp1
integer(ip) :: np1
integer(ip) :: nt
real(wp) :: tsn
real(wp) :: v
real(wp) :: ve
real(wp) :: vo
real(wp) :: w
real(wp) :: we
real(wp) :: wo
real(wp) :: wrfft
real(wp) :: wzvin
real(wp) :: wzwin
real(wp) :: zv
real(wp) :: zw
dimension v(idvw, jdvw, *), w(idvw, jdvw, *), br(mdab, ndab, *), &
          bi(mdab, ndab, *), cr(mdab, ndab, *), ci(mdab, ndab, *), &
          ve(idv, nlon, *), vo(idv, nlon, *), we(idv, nlon, *), &
          wo(idv, nlon, *), wzvin(*), wzwin(*), wrfft(*), &
          zv(imid, nlat, 3), zw(imid, nlat, 3)

type(HFFTpack)      :: hfft
type(SpherepackAux) :: sphere_aux

nlp1 = nlat+1
tsn = 2.0_wp/nlon
fsn = 4.0_wp/nlon
mlat = mod(nlat, 2)
mlon = mod(nlon, 2)
mmax = min(nlat, (nlon+1)/2)
imm1 = imid
if (mlat /= 0) imm1 = imid-1
if (ityp > 2) goto 3  
do 5 k=1, nt 
do 5 i=1, imm1
do 5 j=1, nlon
ve(i, j, k) = tsn*(v(i, j, k)+v(nlp1-i, j, k))
vo(i, j, k) = tsn*(v(i, j, k)-v(nlp1-i, j, k))
we(i, j, k) = tsn*(w(i, j, k)+w(nlp1-i, j, k))
wo(i, j, k) = tsn*(w(i, j, k)-w(nlp1-i, j, k))
5 continue
goto 2
3 do 8 k=1, nt
do 8 i=1, imm1 
do 8 j=1, nlon
ve(i, j, k) = fsn*v(i, j, k)
vo(i, j, k) = fsn*v(i, j, k)
we(i, j, k) = fsn*w(i, j, k)
wo(i, j, k) = fsn*w(i, j, k)
8 continue
2 if (mlat == 0) goto 7
do 6 k=1, nt 
do 6 j=1, nlon
ve(imid, j, k) = tsn*v(imid, j, k)
we(imid, j, k) = tsn*w(imid, j, k)
6 continue
7 do 9 k=1, nt
call hfft%forward(idv, nlon, ve(1, 1, k), idv, wrfft, zv)
call hfft%forward(idv, nlon, we(1, 1, k), idv, wrfft, zv)
9 continue 
ndo1 = nlat
ndo2 = nlat
if (mlat /= 0) ndo1 = nlat-1
if (mlat == 0) ndo2 = nlat-1
if (ityp==2 .or. ityp==5 .or. ityp==8) goto 11 
do 10 k=1, nt
do 10 mp1=1, mmax
do 10 np1=mp1, nlat
br(mp1, np1, k)=0.
bi(mp1, np1, k)=0.
10 continue
11 if (ityp==1 .or. ityp==4 .or. ityp==7) goto 13 
do 12 k=1, nt
do 12 mp1=1, mmax
do 12 np1=mp1, nlat
cr(mp1, np1, k)=0.
ci(mp1, np1, k)=0.
12 continue
13 select case (ityp)
  case (0)
          go to 1
  case (1)
          go to 100
  case (2)
          go to 200
  case (3)
          go to 300
  case (4)
          go to 400
  case (5)
          go to 500
  case (6)
          go to 600
  case (7)
          go to 700
  case (8)
          go to 800
  end select
!
!     case ityp=0 ,  no symmetries
!
1 call sphere_aux%zvin(0, nlat, nlon, 0, zv, iv, wzvin)
!
!     case m=0
!
do 15 k=1, nt
do 15 i=1, imid
do 15 np1=2, ndo2, 2
br(1, np1, k) = br(1, np1, k)+zv(i, np1, iv)*ve(i, 1, k)
cr(1, np1, k) = cr(1, np1, k)-zv(i, np1, iv)*we(i, 1, k)
15 continue
do 16 k=1, nt
do 16 i=1, imm1
do 16 np1=3, ndo1, 2
br(1, np1, k) = br(1, np1, k)+zv(i, np1, iv)*vo(i, 1, k)
cr(1, np1, k) = cr(1, np1, k)-zv(i, np1, iv)*wo(i, 1, k)
16 continue
!
!     case m = 1 through nlat-1
!
if (mmax < 2) return
do 20 mp1=2, mmax
m = mp1-1
mp2 = mp1+1
call sphere_aux%zvin(0, nlat, nlon, m, zv, iv, wzvin)
call sphere_aux%zwin(0, nlat, nlon, m, zw, iw, wzwin)
if (mp1 > ndo1) goto 17
do 23 k=1, nt
do 23 i=1, imm1
do 23 np1=mp1, ndo1, 2
br(mp1, np1, k) = br(mp1, np1, k)+zv(i, np1, iv)*vo(i, 2*mp1-2, k) &
                             +zw(i, np1, iw)*we(i, 2*mp1-1, k)
bi(mp1, np1, k) = bi(mp1, np1, k)+zv(i, np1, iv)*vo(i, 2*mp1-1, k) &
                             -zw(i, np1, iw)*we(i, 2*mp1-2, k)
cr(mp1, np1, k) = cr(mp1, np1, k)-zv(i, np1, iv)*wo(i, 2*mp1-2, k) &
                             +zw(i, np1, iw)*ve(i, 2*mp1-1, k)
ci(mp1, np1, k) = ci(mp1, np1, k)-zv(i, np1, iv)*wo(i, 2*mp1-1, k) &
                             -zw(i, np1, iw)*ve(i, 2*mp1-2, k)
23 continue
if (mlat == 0) goto 17
do 24 k=1, nt
do 24 np1=mp1, ndo1, 2
br(mp1, np1, k) = br(mp1, np1, k)+zw(imid, np1, iw)*we(imid, 2*mp1-1, k)
bi(mp1, np1, k) = bi(mp1, np1, k)-zw(imid, np1, iw)*we(imid, 2*mp1-2, k)
cr(mp1, np1, k) = cr(mp1, np1, k)+zw(imid, np1, iw)*ve(imid, 2*mp1-1, k)
ci(mp1, np1, k) = ci(mp1, np1, k)-zw(imid, np1, iw)*ve(imid, 2*mp1-2, k)
24 continue
17 if (mp2 > ndo2) goto 20
do 21 k=1, nt
do 21 i=1, imm1
do 21 np1=mp2, ndo2, 2
br(mp1, np1, k) = br(mp1, np1, k)+zv(i, np1, iv)*ve(i, 2*mp1-2, k) &
                             +zw(i, np1, iw)*wo(i, 2*mp1-1, k)
bi(mp1, np1, k) = bi(mp1, np1, k)+zv(i, np1, iv)*ve(i, 2*mp1-1, k) &
                             -zw(i, np1, iw)*wo(i, 2*mp1-2, k)
cr(mp1, np1, k) = cr(mp1, np1, k)-zv(i, np1, iv)*we(i, 2*mp1-2, k) &
                             +zw(i, np1, iw)*vo(i, 2*mp1-1, k)
ci(mp1, np1, k) = ci(mp1, np1, k)-zv(i, np1, iv)*we(i, 2*mp1-1, k) &
                             -zw(i, np1, iw)*vo(i, 2*mp1-2, k)
21 continue
if (mlat == 0) goto 20
do 22 k=1, nt
do 22 np1=mp2, ndo2, 2
br(mp1, np1, k) = br(mp1, np1, k)+zv(imid, np1, iv)*ve(imid, 2*mp1-2, k)
bi(mp1, np1, k) = bi(mp1, np1, k)+zv(imid, np1, iv)*ve(imid, 2*mp1-1, k)
cr(mp1, np1, k) = cr(mp1, np1, k)-zv(imid, np1, iv)*we(imid, 2*mp1-2, k)
ci(mp1, np1, k) = ci(mp1, np1, k)-zv(imid, np1, iv)*we(imid, 2*mp1-1, k)
22 continue
20 continue
return
!
!     case ityp=1 ,  no symmetries but cr and ci equal zero
!
100 call sphere_aux%zvin(0, nlat, nlon, 0, zv, iv, wzvin)
!
!     case m=0
!
do 115 k=1, nt
do 115 i=1, imid
do 115 np1=2, ndo2, 2
br(1, np1, k) = br(1, np1, k)+zv(i, np1, iv)*ve(i, 1, k)
115 continue
do 116 k=1, nt
do 116 i=1, imm1
do 116 np1=3, ndo1, 2
br(1, np1, k) = br(1, np1, k)+zv(i, np1, iv)*vo(i, 1, k)
116 continue
!
!     case m = 1 through nlat-1
!
if (mmax < 2) return
do 120 mp1=2, mmax
m = mp1-1
mp2 = mp1+1
call sphere_aux%zvin(0, nlat, nlon, m, zv, iv, wzvin)
call sphere_aux%zwin(0, nlat, nlon, m, zw, iw, wzwin)
if (mp1 > ndo1) goto 117
do 123 k=1, nt
do 123 i=1, imm1
do 123 np1=mp1, ndo1, 2
br(mp1, np1, k) = br(mp1, np1, k)+zv(i, np1, iv)*vo(i, 2*mp1-2, k) &
                             +zw(i, np1, iw)*we(i, 2*mp1-1, k)
bi(mp1, np1, k) = bi(mp1, np1, k)+zv(i, np1, iv)*vo(i, 2*mp1-1, k) &
                             -zw(i, np1, iw)*we(i, 2*mp1-2, k)
123 continue
if (mlat == 0) goto 117
do 124 k=1, nt
do 124 np1=mp1, ndo1, 2
br(mp1, np1, k) = br(mp1, np1, k)+zw(imid, np1, iw)*we(imid, 2*mp1-1, k)
bi(mp1, np1, k) = bi(mp1, np1, k)-zw(imid, np1, iw)*we(imid, 2*mp1-2, k)
124 continue
117 if (mp2 > ndo2) goto 120
do 121 k=1, nt
do 121 i=1, imm1
do 121 np1=mp2, ndo2, 2
br(mp1, np1, k) = br(mp1, np1, k)+zv(i, np1, iv)*ve(i, 2*mp1-2, k) &
                             +zw(i, np1, iw)*wo(i, 2*mp1-1, k)
bi(mp1, np1, k) = bi(mp1, np1, k)+zv(i, np1, iv)*ve(i, 2*mp1-1, k) &
                             -zw(i, np1, iw)*wo(i, 2*mp1-2, k)
121 continue
if (mlat == 0) goto 120
do 122 k=1, nt
do 122 np1=mp2, ndo2, 2
br(mp1, np1, k) = br(mp1, np1, k)+zv(imid, np1, iv)*ve(imid, 2*mp1-2, k)
bi(mp1, np1, k) = bi(mp1, np1, k)+zv(imid, np1, iv)*ve(imid, 2*mp1-1, k)
122 continue
120 continue
return
!
!     case ityp=2 ,  no symmetries but br and bi equal zero   
!
200 call sphere_aux%zvin(0, nlat, nlon, 0, zv, iv, wzvin)
!
!     case m=0
!
do 215 k=1, nt
do 215 i=1, imid
do 215 np1=2, ndo2, 2
cr(1, np1, k) = cr(1, np1, k)-zv(i, np1, iv)*we(i, 1, k)
215 continue
do 216 k=1, nt
do 216 i=1, imm1
do 216 np1=3, ndo1, 2
cr(1, np1, k) = cr(1, np1, k)-zv(i, np1, iv)*wo(i, 1, k)
216 continue
!
!     case m = 1 through nlat-1
!
if (mmax < 2) return
do 220 mp1=2, mmax
m = mp1-1
mp2 = mp1+1
call sphere_aux%zvin(0, nlat, nlon, m, zv, iv, wzvin)
call sphere_aux%zwin(0, nlat, nlon, m, zw, iw, wzwin)
if (mp1 > ndo1) goto 217
do 223 k=1, nt
do 223 i=1, imm1
do 223 np1=mp1, ndo1, 2
cr(mp1, np1, k) = cr(mp1, np1, k)-zv(i, np1, iv)*wo(i, 2*mp1-2, k) &
                             +zw(i, np1, iw)*ve(i, 2*mp1-1, k)
ci(mp1, np1, k) = ci(mp1, np1, k)-zv(i, np1, iv)*wo(i, 2*mp1-1, k) &
                             -zw(i, np1, iw)*ve(i, 2*mp1-2, k)
223 continue
if (mlat == 0) goto 217
do 224 k=1, nt
do 224 np1=mp1, ndo1, 2
cr(mp1, np1, k) = cr(mp1, np1, k)+zw(imid, np1, iw)*ve(imid, 2*mp1-1, k)
ci(mp1, np1, k) = ci(mp1, np1, k)-zw(imid, np1, iw)*ve(imid, 2*mp1-2, k)
224 continue
217 if (mp2 > ndo2) goto 220
do 221 k=1, nt
do 221 i=1, imm1
do 221 np1=mp2, ndo2, 2
cr(mp1, np1, k) = cr(mp1, np1, k)-zv(i, np1, iv)*we(i, 2*mp1-2, k) &
                             +zw(i, np1, iw)*vo(i, 2*mp1-1, k)
ci(mp1, np1, k) = ci(mp1, np1, k)-zv(i, np1, iv)*we(i, 2*mp1-1, k) &
                             -zw(i, np1, iw)*vo(i, 2*mp1-2, k)
221 continue
if (mlat == 0) goto 220
do 222 k=1, nt
do 222 np1=mp2, ndo2, 2
cr(mp1, np1, k) = cr(mp1, np1, k)-zv(imid, np1, iv)*we(imid, 2*mp1-2, k)
ci(mp1, np1, k) = ci(mp1, np1, k)-zv(imid, np1, iv)*we(imid, 2*mp1-1, k)
222 continue
220 continue
return
!
!     case ityp=3 ,  v even , w odd
!
300 call sphere_aux%zvin(0, nlat, nlon, 0, zv, iv, wzvin)
!
!     case m=0
!
do 315 k=1, nt
do 315 i=1, imid
do 315 np1=2, ndo2, 2
br(1, np1, k) = br(1, np1, k)+zv(i, np1, iv)*ve(i, 1, k)
315 continue
do 316 k=1, nt
do 316 i=1, imm1
do 316 np1=3, ndo1, 2
cr(1, np1, k) = cr(1, np1, k)-zv(i, np1, iv)*wo(i, 1, k)
316 continue
!
!     case m = 1 through nlat-1
!
if (mmax < 2) return
do 320 mp1=2, mmax
m = mp1-1
mp2 = mp1+1
call sphere_aux%zvin(0, nlat, nlon, m, zv, iv, wzvin)
call sphere_aux%zwin(0, nlat, nlon, m, zw, iw, wzwin)
if (mp1 > ndo1) goto 317
do 323 k=1, nt
do 323 i=1, imm1
do 323 np1=mp1, ndo1, 2
cr(mp1, np1, k) = cr(mp1, np1, k)-zv(i, np1, iv)*wo(i, 2*mp1-2, k) &
                             +zw(i, np1, iw)*ve(i, 2*mp1-1, k)
ci(mp1, np1, k) = ci(mp1, np1, k)-zv(i, np1, iv)*wo(i, 2*mp1-1, k) &
                             -zw(i, np1, iw)*ve(i, 2*mp1-2, k)
323 continue
if (mlat == 0) goto 317
do 324 k=1, nt
do 324 np1=mp1, ndo1, 2
cr(mp1, np1, k) = cr(mp1, np1, k)+zw(imid, np1, iw)*ve(imid, 2*mp1-1, k)
ci(mp1, np1, k) = ci(mp1, np1, k)-zw(imid, np1, iw)*ve(imid, 2*mp1-2, k)
324 continue
317 if (mp2 > ndo2) goto 320
do 321 k=1, nt
do 321 i=1, imm1
do 321 np1=mp2, ndo2, 2
br(mp1, np1, k) = br(mp1, np1, k)+zv(i, np1, iv)*ve(i, 2*mp1-2, k) &
                             +zw(i, np1, iw)*wo(i, 2*mp1-1, k)
bi(mp1, np1, k) = bi(mp1, np1, k)+zv(i, np1, iv)*ve(i, 2*mp1-1, k) &
                             -zw(i, np1, iw)*wo(i, 2*mp1-2, k)
321 continue
if (mlat == 0) goto 320
do 322 k=1, nt
do 322 np1=mp2, ndo2, 2
br(mp1, np1, k) = br(mp1, np1, k)+zv(imid, np1, iv)*ve(imid, 2*mp1-2, k)
bi(mp1, np1, k) = bi(mp1, np1, k)+zv(imid, np1, iv)*ve(imid, 2*mp1-1, k)
322 continue
320 continue
return
!
!     case ityp=4 ,  v even, w odd, and cr and ci equal 0. 
!
400 call sphere_aux%zvin(1, nlat, nlon, 0, zv, iv, wzvin)
!
!     case m=0
!
do 415 k=1, nt
do 415 i=1, imid
do 415 np1=2, ndo2, 2
br(1, np1, k) = br(1, np1, k)+zv(i, np1, iv)*ve(i, 1, k)
415 continue
!
!     case m = 1 through nlat-1
!
if (mmax < 2) return
do 420 mp1=2, mmax
m = mp1-1
mp2 = mp1+1
call sphere_aux%zvin(1, nlat, nlon, m, zv, iv, wzvin)
call sphere_aux%zwin(1, nlat, nlon, m, zw, iw, wzwin)
if (mp2 > ndo2) goto 420
do 421 k=1, nt
do 421 i=1, imm1
do 421 np1=mp2, ndo2, 2
br(mp1, np1, k) = br(mp1, np1, k)+zv(i, np1, iv)*ve(i, 2*mp1-2, k) &
                             +zw(i, np1, iw)*wo(i, 2*mp1-1, k)
bi(mp1, np1, k) = bi(mp1, np1, k)+zv(i, np1, iv)*ve(i, 2*mp1-1, k) &
                             -zw(i, np1, iw)*wo(i, 2*mp1-2, k)
421 continue
if (mlat == 0) goto 420
do 422 k=1, nt
do 422 np1=mp2, ndo2, 2
br(mp1, np1, k) = br(mp1, np1, k)+zv(imid, np1, iv)*ve(imid, 2*mp1-2, k)
bi(mp1, np1, k) = bi(mp1, np1, k)+zv(imid, np1, iv)*ve(imid, 2*mp1-1, k)
422 continue
420 continue
return
!
!     case ityp=5   v even, w odd, and br and bi equal zero
!
500 call sphere_aux%zvin(2, nlat, nlon, 0, zv, iv, wzvin)
!
!     case m=0
!
do 516 k=1, nt
do 516 i=1, imm1
do 516 np1=3, ndo1, 2
cr(1, np1, k) = cr(1, np1, k)-zv(i, np1, iv)*wo(i, 1, k)
516 continue
!
!     case m = 1 through nlat-1
!
if (mmax < 2) return
do 520 mp1=2, mmax
m = mp1-1
mp2 = mp1+1
call sphere_aux%zvin(2, nlat, nlon, m, zv, iv, wzvin)
call sphere_aux%zwin(2, nlat, nlon, m, zw, iw, wzwin)
if (mp1 > ndo1) goto 520
do 523 k=1, nt
do 523 i=1, imm1
do 523 np1=mp1, ndo1, 2
cr(mp1, np1, k) = cr(mp1, np1, k)-zv(i, np1, iv)*wo(i, 2*mp1-2, k) &
                             +zw(i, np1, iw)*ve(i, 2*mp1-1, k)
ci(mp1, np1, k) = ci(mp1, np1, k)-zv(i, np1, iv)*wo(i, 2*mp1-1, k) &
                             -zw(i, np1, iw)*ve(i, 2*mp1-2, k)
523 continue
if (mlat == 0) goto 520
do 524 k=1, nt
do 524 np1=mp1, ndo1, 2
cr(mp1, np1, k) = cr(mp1, np1, k)+zw(imid, np1, iw)*ve(imid, 2*mp1-1, k)
ci(mp1, np1, k) = ci(mp1, np1, k)-zw(imid, np1, iw)*ve(imid, 2*mp1-2, k)
524 continue
520 continue
return
!
!     case ityp=6 ,  v odd , w even
!
600 call sphere_aux%zvin(0, nlat, nlon, 0, zv, iv, wzvin)
!
!     case m=0
!
do 615 k=1, nt
do 615 i=1, imid
do 615 np1=2, ndo2, 2
cr(1, np1, k) = cr(1, np1, k)-zv(i, np1, iv)*we(i, 1, k)
615 continue
do 616 k=1, nt
do 616 i=1, imm1
do 616 np1=3, ndo1, 2
br(1, np1, k) = br(1, np1, k)+zv(i, np1, iv)*vo(i, 1, k)
616 continue
!
!     case m = 1 through nlat-1
!
if (mmax < 2) return
do 620 mp1=2, mmax
m = mp1-1
mp2 = mp1+1
call sphere_aux%zvin(0, nlat, nlon, m, zv, iv, wzvin)
call sphere_aux%zwin(0, nlat, nlon, m, zw, iw, wzwin)
if (mp1 > ndo1) goto 617
do 623 k=1, nt
do 623 i=1, imm1
do 623 np1=mp1, ndo1, 2
br(mp1, np1, k) = br(mp1, np1, k)+zv(i, np1, iv)*vo(i, 2*mp1-2, k) &
                             +zw(i, np1, iw)*we(i, 2*mp1-1, k)
bi(mp1, np1, k) = bi(mp1, np1, k)+zv(i, np1, iv)*vo(i, 2*mp1-1, k) &
                             -zw(i, np1, iw)*we(i, 2*mp1-2, k)
623 continue
if (mlat == 0) goto 617
do 624 k=1, nt
do 624 np1=mp1, ndo1, 2
br(mp1, np1, k) = br(mp1, np1, k)+zw(imid, np1, iw)*we(imid, 2*mp1-1, k)
bi(mp1, np1, k) = bi(mp1, np1, k)-zw(imid, np1, iw)*we(imid, 2*mp1-2, k)
624 continue
617 if (mp2 > ndo2) goto 620
do 621 k=1, nt
do 621 i=1, imm1
do 621 np1=mp2, ndo2, 2
cr(mp1, np1, k) = cr(mp1, np1, k)-zv(i, np1, iv)*we(i, 2*mp1-2, k) &
                             +zw(i, np1, iw)*vo(i, 2*mp1-1, k)
ci(mp1, np1, k) = ci(mp1, np1, k)-zv(i, np1, iv)*we(i, 2*mp1-1, k) &
                             -zw(i, np1, iw)*vo(i, 2*mp1-2, k)
621 continue
if (mlat == 0) goto 620
do 622 k=1, nt
do 622 np1=mp2, ndo2, 2
cr(mp1, np1, k) = cr(mp1, np1, k)-zv(imid, np1, iv)*we(imid, 2*mp1-2, k)
ci(mp1, np1, k) = ci(mp1, np1, k)-zv(imid, np1, iv)*we(imid, 2*mp1-1, k)
622 continue
620 continue
return
!
!     case ityp=7   v odd, w even, and cr and ci equal zero
!
700 call sphere_aux%zvin(2, nlat, nlon, 0, zv, iv, wzvin)
!
!     case m=0
!
do 716 k=1, nt
do 716 i=1, imm1
do 716 np1=3, ndo1, 2
br(1, np1, k) = br(1, np1, k)+zv(i, np1, iv)*vo(i, 1, k)
716 continue
!
!     case m = 1 through nlat-1
!
if (mmax < 2) return
do 720 mp1=2, mmax
m = mp1-1
mp2 = mp1+1
call sphere_aux%zvin(2, nlat, nlon, m, zv, iv, wzvin)
call sphere_aux%zwin(2, nlat, nlon, m, zw, iw, wzwin)
if (mp1 > ndo1) goto 720
do k=1, nt
do i=1, imm1
do np1=mp1, ndo1, 2
br(mp1, np1, k) = br(mp1, np1, k)+zv(i, np1, iv)*vo(i, 2*mp1-2, k) &
                             +zw(i, np1, iw)*we(i, 2*mp1-1, k)
bi(mp1, np1, k) = bi(mp1, np1, k)+zv(i, np1, iv)*vo(i, 2*mp1-1, k) &
                             -zw(i, np1, iw)*we(i, 2*mp1-2, k)
end do
end do
end do
if (mlat == 0) goto 720
do k=1, nt
do np1=mp1, ndo1, 2
br(mp1, np1, k) = br(mp1, np1, k)+zw(imid, np1, iw)*we(imid, 2*mp1-1, k)
bi(mp1, np1, k) = bi(mp1, np1, k)-zw(imid, np1, iw)*we(imid, 2*mp1-2, k)
end do
end do
720 continue
return
!
!     case ityp=8   v odd, w even, and both br and bi equal zero
!
800 call sphere_aux%zvin(1, nlat, nlon, 0, zv, iv, wzvin)
!
!     case m=0
!
do 815 k=1, nt
do 815 i=1, imid
do 815 np1=2, ndo2, 2
cr(1, np1, k) = cr(1, np1, k)-zv(i, np1, iv)*we(i, 1, k)
815 continue
!
!     case m = 1 through nlat-1
!
if (mmax < 2) return
do 820 mp1=2, mmax
m = mp1-1
mp2 = mp1+1
call sphere_aux%zvin(1, nlat, nlon, m, zv, iv, wzvin)
call sphere_aux%zwin(1, nlat, nlon, m, zw, iw, wzwin)
if (mp2 > ndo2) goto 820
do k=1, nt
do i=1, imm1
do np1=mp2, ndo2, 2
cr(mp1, np1, k) = cr(mp1, np1, k)-zv(i, np1, iv)*we(i, 2*mp1-2, k) &
                             +zw(i, np1, iw)*vo(i, 2*mp1-1, k)
ci(mp1, np1, k) = ci(mp1, np1, k)-zv(i, np1, iv)*we(i, 2*mp1-1, k) &
                             -zw(i, np1, iw)*vo(i, 2*mp1-2, k)
end do
end do
end do
if (mlat == 0) goto 820
do k=1, nt
do np1=mp2, ndo2, 2
cr(mp1, np1, k) = cr(mp1, np1, k)-zv(imid, np1, iv)*we(imid, 2*mp1-2, k)
ci(mp1, np1, k) = ci(mp1, np1, k)-zv(imid, np1, iv)*we(imid, 2*mp1-1, k)
end do
end do
820 continue

end subroutine vhaec1

end subroutine vhaec


subroutine vhaeci(nlat, nlon, wvhaec, lvhaec, dwork, ldwork, ierror)

integer(ip) :: ierror
integer(ip) :: imid
integer(ip) :: iw1
integer(ip) :: iw2
integer(ip) :: labc
integer(ip) :: ldwork
integer(ip) :: lvhaec
integer(ip) :: lwzvin
integer(ip) :: lzz1
integer(ip) :: mmax
integer(ip) :: nlat
integer(ip) :: nlon
real(wp) :: wvhaec(lvhaec)
real(wp) :: dwork(ldwork)

type(HFFTpack)      :: hfft
type(SpherepackAux) :: sphere_aux

imid = (nlat+1)/2
lzz1 = 2*nlat*imid
mmax = min(nlat, (nlon+1)/2)
labc = 3*(max(mmax-2, 0)*(2*nlat-mmax-1))/2


ierror = 1
if (nlat < 3) return
ierror = 2
if (nlon < 1) return
ierror = 3

if (lvhaec < 2*(lzz1+labc)+nlon+15) return
ierror = 4
if (ldwork < 2*nlat+2) return
ierror = 0

call sphere_aux%zvinit(nlat, nlon, wvhaec, dwork)

!
!==> Set workspace pointers
!
lwzvin = lzz1+labc
iw1 = lwzvin+1
iw2 = iw1+lwzvin

call sphere_aux%zwinit(nlat, nlon, wvhaec(iw1), dwork)

call hfft%initialize(nlon, wvhaec(iw2))

end subroutine vhaeci

end module module_vhaec
