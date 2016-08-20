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
! ... file vhagc.f
!
!     this file contains code and documentation for subroutines
!     vhagc and vhagci
!
! ... files which must be loaded with vhagc.f
!
!     type_SpherepackAux.f, type_HFFTpack.f, compute_gaussian_latitudes_and_weights.f
!
!                                                                              
!     subroutine vhagc(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, 
!    +                 mdab, ndab, wvhagc, lvhagc, work, lwork, ierror)
!
!     subroutine vhagc performs the vector spherical harmonic analysis
!     on the vector field (v, w) and stores the result in the arrays
!     br, bi, cr, and ci. v(i, j) and w(i, j) are the colatitudinal
!     (measured from the north pole) and east longitudinal components
!     respectively, located at the gaussian colatitude point theta(i)
!     and longitude phi(j) = (j-1)*2*pi/nlon. the spectral
!     representation of (v, w) is given at output parameters v, w in 
!     subroutine vhsec.  
!
!     input parameters
!
!     nlat   the number of points in the gaussian colatitude grid on the
!            full sphere. these lie in the interval (0, pi) and are computed
!            in radians in theta(1) <...< theta(nlat) by subroutine compute_gaussian_latitudes_and_weights.
!            if nlat is odd the equator will be included as the grid point
!            theta((nlat+1)/2).  if nlat is even the equator will be
!            excluded as a grid point and will lie half way between
!            theta(nlat/2) and theta(nlat/2+1). nlat must be at least 3.
!            note: on the half sphere, the number of grid points in the
!            colatitudinal direction is nlat/2 if nlat is even or
!            (nlat+1)/2 if nlat is odd.
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
!     nt     the number of analyses.  in the program that calls vhagc, 
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
!            the program that calls vhagc. if ityp .le. 2 then idvw
!            must be at least nlat.  if ityp .gt. 2 and nlat is
!            even then idvw must be at least nlat/2. if ityp .gt. 2
!            and nlat is odd then idvw must be at least (nlat+1)/2.
!
!     jdvw   the second dimension of the arrays v, w as it appears in
!            the program that calls vhagc. jdvw must be at least nlon.
!
!     mdab   the first dimension of the arrays br, bi, cr, and ci as it
!            appears in the program that calls vhagc. mdab must be at
!            least min(nlat, nlon/2) if nlon is even or at least
!            min(nlat, (nlon+1)/2) if nlon is odd.
!
!     ndab   the second dimension of the arrays br, bi, cr, and ci as it
!            appears in the program that calls vhagc. ndab must be at
!            least nlat.
!
!     wvhagc an array which must be initialized by subroutine vhagci.
!            once initialized, wvhagc can be used repeatedly by vhagc
!            as long as nlon and nlat remain unchanged.  wvhagc must
!            not be altered between calls of vhagc.
!
!     lvhagc the dimension of the array wvhagc as it appears in the
!            program that calls vhagc. define
!
!               l1 = min(nlat, nlon/2) if nlon is even or
!               l1 = min(nlat, (nlon+1)/2) if nlon is odd
!
!            and
!
!               l2 = nlat/2        if nlat is even or
!               l2 = (nlat+1)/2    if nlat is odd
!
!            then lvhagc must be at least
!
!               4*nlat*l2+3*max(l1-2, 0)*(2*nlat-l1-1)+nlon+l2+15
!
!
!     work   a work array that does not have to be saved.
!
!     lwork  the dimension of the array work as it appears in the
!            program that calls vhagc. define
!
!               l2 = nlat/2        if nlat is even or
!               l2 = (nlat+1)/2    if nlat is odd
!
!            if ityp .le. 2 then lwork must be at least
!
!               2*nlat*(2*nlon*nt+3*l2)
!
!            if ityp .gt. 2 then lwork must be at least
!
!               2*l2*(2*nlon*nt+3*nlat)
!
!
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
!            = 9  error in the specification of lvhagc
!            = 10 error in the specification of lwork
!
! ****************************************************************
!
!     subroutine vhagci(nlat, nlon, wvhagc, lvhagc, dwork, ldwork, ierror)
!
!     subroutine vhagci initializes the array wvhagc which can then be
!     used repeatedly by subroutine vhagc until nlat or nlon is changed.
!
!     input parameters
!
!     nlat   the number of points in the gaussian colatitude grid on the
!            full sphere. these lie in the interval (0, pi) and are computed
!            in radians in theta(1) <...< theta(nlat) by subroutine compute_gaussian_latitudes_and_weights.
!            if nlat is odd the equator will be included as the grid point
!            theta((nlat+1)/2).  if nlat is even the equator will be
!            excluded as a grid point and will lie half way between
!            theta(nlat/2) and theta(nlat/2+1). nlat must be at least 3.
!            note: on the half sphere, the number of grid points in the
!            colatitudinal direction is nlat/2 if nlat is even or
!            (nlat+1)/2 if nlat is odd.
!
!     nlon   the number of distinct londitude points.  nlon determines
!            the grid increment in longitude as 2*pi/nlon. for example
!            nlon = 72 for a five degree grid. nlon must be greater
!            than zero. the axisymmetric case corresponds to nlon=1.
!            the efficiency of the computation is improved when nlon
!            is a product of small prime numbers.
!
!     lvhagc the dimension of the array wvhagc as it appears in the
!            program that calls vhagci.  define
!
!               l1 = min(nlat, nlon/2) if nlon is even or
!               l1 = min(nlat, (nlon+1)/2) if nlon is odd
!
!            and
!
!               l2 = nlat/2        if nlat is even or
!               l2 = (nlat+1)/2    if nlat is odd
!
!            then lvhagc must be at least
!
!               4*nlat*l2+3*max(l1-2, 0)*(2*nlat-l1-1)+nlon+l2+15
!
!
!     dwork  a real work array that does not need to be saved
!
!     ldwork the dimension of the array dwork as it appears in the
!            program that calls vhagci. ldwork must be at least
!
!               2*nlat*(nlat+1)+1
!
!
!     **************************************************************
!
!     output parameters
!
!     wvhagc an array which is initialized for use by subroutine vhagc.
!            once initialized, wvhagc can be used repeatedly by vhagc
!            as long as nlat and nlon remain unchanged.  wvhagc must not
!            be altered between calls of vhagc.
!
!
!     ierror = 0  no errors
!            = 1  error in the specification of nlat
!            = 2  error in the specification of nlon
!            = 3  error in the specification of lvhagc
!            = 4  error in the specification of lwork
!
module module_vhagc

    use spherepack_precision, only: &
        wp, & ! working precision
        ip ! integer precision

    use type_HFFTpack, only: &
        HFFTpack

    use type_SpherepackAux, only: &
        SpherepackAux

    use gaussian_latitudes_and_weights_routines, only: &
        compute_gaussian_latitudes_and_weights

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    public :: vhagc
    public :: vhagci

contains

subroutine vhagc(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, &
           mdab, ndab, wvhagc, lvhagc, work, lwork, ierror)

real (wp) :: bi
real (wp) :: br
real (wp) :: ci
real (wp) :: cr
integer (ip) :: idv
integer (ip) :: idvw
integer (ip) :: ierror
integer (ip) :: imid
integer (ip) :: ist
integer (ip) :: ityp
integer (ip) :: iw1
integer (ip) :: iw2
integer (ip) :: iw3
integer (ip) :: iw4
integer (ip) :: iw5
integer (ip) :: jdvw
integer (ip) :: jw1
integer (ip) :: jw2
integer (ip) :: jw3
integer (ip) :: labc
integer (ip) :: lnl
integer (ip) :: lvhagc
integer (ip) :: lwork
integer (ip) :: lwzvin
integer (ip) :: lzz1
integer (ip) :: mdab
integer (ip) :: mmax
integer (ip) :: ndab
integer (ip) :: nlat
integer (ip) :: nlon
integer (ip) :: nt
real (wp) :: v
real (wp) :: w
real (wp) :: work(lwork)
real (wp) :: wvhagc(lvhagc)
dimension v(idvw, jdvw, *), w(idvw, jdvw, *), br(mdab, ndab, *), &
          bi(mdab, ndab, *), cr(mdab, ndab, *), ci(mdab, ndab, *)

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
if (lvhagc < 2*(lzz1+labc)+nlon+imid+15) return
ierror = 10
if (ityp<=2 .and. lwork <nlat*(4*nlon*nt+6*imid)) return
if (ityp>2 .and. lwork <imid*(4*nlon*nt+6*nlat)) return
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
jw1 = (nlat+1)/2+1
jw2 = jw1+lwzvin
jw3 = jw2+lwzvin
call vhagc1(nlat, nlon, ityp, nt, imid, idvw, jdvw, v, w, mdab, ndab, &
br, bi, cr, ci, idv, work, work(iw1), work(iw2), work(iw3), &
work(iw4), work(iw5), wvhagc, wvhagc(jw1), wvhagc(jw2), wvhagc(jw3))

contains

subroutine vhagc1(nlat, nlon, ityp, nt, imid, idvw, jdvw, v, w, mdab, &
ndab, br, bi, cr, ci, idv, ve, vo, we, wo, vb, wb, wts, wvbin, wwbin, wrfft)

real (wp) :: bi
real (wp) :: br
real (wp) :: ci
real (wp) :: cr
real (wp) :: fsn
integer (ip) :: i
integer (ip) :: idv
integer (ip) :: idvw
integer (ip) :: imid
integer (ip) :: imm1
integer (ip) :: ityp
integer (ip) :: itypp
integer (ip) :: iv
integer (ip) :: iw
integer (ip) :: j
integer (ip) :: jdvw
integer (ip) :: k
integer (ip) :: m
integer (ip) :: mdab
integer (ip) :: mlat
integer (ip) :: mlon
integer (ip) :: mmax
integer (ip) :: mp1
integer (ip) :: mp2
integer (ip) :: ndab
integer (ip) :: ndo1
integer (ip) :: ndo2
integer (ip) :: nlat
integer (ip) :: nlon
integer (ip) :: nlp1
integer (ip) :: np1
integer (ip) :: nt
real (wp) :: tsn
real (wp) :: tv
real (wp) :: tve1
real (wp) :: tve2
real (wp) :: tvo1
real (wp) :: tvo2
real (wp) :: tw
real (wp) :: twe1
real (wp) :: twe2
real (wp) :: two1
real (wp) :: two2
real (wp) :: v
real (wp) :: vb
real (wp) :: ve
real (wp) :: vo
real (wp) :: w
real (wp) :: wb
real (wp) :: we
real (wp) :: wo
real (wp) :: wrfft
real (wp) :: wts
real (wp) :: wvbin
real (wp) :: wwbin
dimension v(idvw, jdvw, *), w(idvw, jdvw, *), br(mdab, ndab, *), &
          bi(mdab, ndab, *), cr(mdab, ndab, *), ci(mdab, ndab, *), &
          ve(idv, nlon, *), vo(idv, nlon, *), we(idv, nlon, *), &
          wo(idv, nlon, *), wts(*), wvbin(*), wwbin(*), wrfft(*), &
          vb(imid, nlat, 3), wb(imid, nlat, 3)

type (HFFTpack)      :: hfft
type (SpherepackAux) :: sphere_aux

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
call hfft%forward(idv, nlon, ve(1, 1, k), idv, wrfft, vb)
call hfft%forward(idv, nlon, we(1, 1, k), idv, wrfft, vb)
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
13 itypp = ityp+1
goto (1, 100, 200, 300, 400, 500, 600, 700, 800), itypp
!
!     case ityp=0 ,  no symmetries
!
1 call sphere_aux%vbin(0, nlat, nlon, 0, vb, iv, wvbin)
!
!     case m=0
!
do 15 k=1, nt
do 1015 i=1, imid
tv = ve(i, 1, k)*wts(i)
tw = we(i, 1, k)*wts(i)
do 10015 np1=2, ndo2, 2
br(1, np1, k) = br(1, np1, k)+vb(i, np1, iv)*tv
cr(1, np1, k) = cr(1, np1, k)-vb(i, np1, iv)*tw
10015 continue
1015 continue
15 continue
do 16 k=1, nt
do 1016 i=1, imm1
tv = vo(i, 1, k)*wts(i)
tw = wo(i, 1, k)*wts(i)
do 10016 np1=3, ndo1, 2
br(1, np1, k) = br(1, np1, k)+vb(i, np1, iv)*tv
cr(1, np1, k) = cr(1, np1, k)-vb(i, np1, iv)*tw
10016 continue
1016 continue
16 continue
!
!     case m = 1 through nlat-1
!
if (mmax < 2) return
do 20 mp1=2, mmax
m = mp1-1
mp2 = mp1+1
call sphere_aux%vbin(0, nlat, nlon, m, vb, iv, wvbin)
call sphere_aux%wbin(0, nlat, nlon, m, wb, iw, wwbin)
if (mp1 > ndo1) goto 17
do 23 k=1, nt
do 1023 i=1, imm1
!
!     set temps to optimize quadrature
!
tvo1 = vo(i, 2*mp1-1, k)*wts(i)
tvo2 = vo(i, 2*mp1-2, k)*wts(i)
tve1 = ve(i, 2*mp1-1, k)*wts(i)
tve2 = ve(i, 2*mp1-2, k)*wts(i)
two1 = wo(i, 2*mp1-1, k)*wts(i)
two2 = wo(i, 2*mp1-2, k)*wts(i)
twe1 = we(i, 2*mp1-1, k)*wts(i)
twe2 = we(i, 2*mp1-2, k)*wts(i)
do 10023 np1=mp1, ndo1, 2
br(mp1, np1, k) = br(mp1, np1, k)+vb(i, np1, iv)*tvo2 &
                             +wb(i, np1, iw)*twe1
bi(mp1, np1, k) = bi(mp1, np1, k)+vb(i, np1, iv)*tvo1 &
                             -wb(i, np1, iw)*twe2
cr(mp1, np1, k) = cr(mp1, np1, k)-vb(i, np1, iv)*two2 &
                             +wb(i, np1, iw)*tve1
ci(mp1, np1, k) = ci(mp1, np1, k)-vb(i, np1, iv)*two1 &
                             -wb(i, np1, iw)*tve2
10023 continue
1023 continue
23 continue

if (mlat == 0) goto 17
i = imid
do 24 k=1, nt
do 1024 np1=mp1, ndo1, 2
br(mp1, np1, k)=br(mp1, np1, k)+wb(i, np1, iw)*we(i, 2*mp1-1, k)*wts(i)
bi(mp1, np1, k)=bi(mp1, np1, k)-wb(i, np1, iw)*we(i, 2*mp1-2, k)*wts(i)
cr(mp1, np1, k)=cr(mp1, np1, k)+wb(i, np1, iw)*ve(i, 2*mp1-1, k)*wts(i)
ci(mp1, np1, k)=ci(mp1, np1, k)-wb(i, np1, iw)*ve(i, 2*mp1-2, k)*wts(i)
1024 continue
24 continue
17 if (mp2 > ndo2) goto 20
do 21 k=1, nt
do 1021 i=1, imm1
tvo1 = vo(i, 2*mp1-1, k)*wts(i)
tvo2 = vo(i, 2*mp1-2, k)*wts(i)
tve1 = ve(i, 2*mp1-1, k)*wts(i)
tve2 = ve(i, 2*mp1-2, k)*wts(i)
two1 = wo(i, 2*mp1-1, k)*wts(i)
two2 = wo(i, 2*mp1-2, k)*wts(i)
twe1 = we(i, 2*mp1-1, k)*wts(i)
twe2 = we(i, 2*mp1-2, k)*wts(i)
do 10021 np1=mp2, ndo2, 2
br(mp1, np1, k) = br(mp1, np1, k)+vb(i, np1, iv)*tve2 &
                             +wb(i, np1, iw)*two1
bi(mp1, np1, k) = bi(mp1, np1, k)+vb(i, np1, iv)*tve1 &
                             -wb(i, np1, iw)*two2
cr(mp1, np1, k) = cr(mp1, np1, k)-vb(i, np1, iv)*twe2 &
                             +wb(i, np1, iw)*tvo1
ci(mp1, np1, k) = ci(mp1, np1, k)-vb(i, np1, iv)*twe1 &
                             -wb(i, np1, iw)*tvo2
10021 continue
1021 continue
21 continue

if (mlat == 0) goto 20
i = imid
do 22 k=1, nt
do 1022 np1=mp2, ndo2, 2
br(mp1, np1, k)=br(mp1, np1, k)+vb(i, np1, iv)*ve(i, 2*mp1-2, k)*wts(i)
bi(mp1, np1, k)=bi(mp1, np1, k)+vb(i, np1, iv)*ve(i, 2*mp1-1, k)*wts(i)
cr(mp1, np1, k)=cr(mp1, np1, k)-vb(i, np1, iv)*we(i, 2*mp1-2, k)*wts(i)
ci(mp1, np1, k)=ci(mp1, np1, k)-vb(i, np1, iv)*we(i, 2*mp1-1, k)*wts(i)
1022 continue
22 continue
20 continue
return
!
!     case ityp=1 ,  no symmetries but cr and ci equal zero
!
100 call sphere_aux%vbin(0, nlat, nlon, 0, vb, iv, wvbin)
!
!     case m=0
!
do 115 k=1, nt
do 115 i=1, imid
tv = ve(i, 1, k)*wts(i)
do 115 np1=2, ndo2, 2
br(1, np1, k) = br(1, np1, k)+vb(i, np1, iv)*tv
115 continue
do 116 k=1, nt
do 116 i=1, imm1
tv = vo(i, 1, k)*wts(i)
do 116 np1=3, ndo1, 2
br(1, np1, k) = br(1, np1, k)+vb(i, np1, iv)*tv
116 continue
!
!     case m = 1 through nlat-1
!
if (mmax < 2) return
do 120 mp1=2, mmax
m = mp1-1
mp2 = mp1+1
call sphere_aux%vbin(0, nlat, nlon, m, vb, iv, wvbin)
call sphere_aux%wbin(0, nlat, nlon, m, wb, iw, wwbin)
if (mp1 > ndo1) goto 117
do 123 k=1, nt
do 123 i=1, imm1
tvo1 = vo(i, 2*mp1-1, k)*wts(i)
tvo2 = vo(i, 2*mp1-2, k)*wts(i)
twe1 = we(i, 2*mp1-1, k)*wts(i)
twe2 = we(i, 2*mp1-2, k)*wts(i)
do 123 np1=mp1, ndo1, 2
br(mp1, np1, k) = br(mp1, np1, k)+vb(i, np1, iv)*tvo2 &
                             +wb(i, np1, iw)*twe1
bi(mp1, np1, k) = bi(mp1, np1, k)+vb(i, np1, iv)*tvo1 &
                             -wb(i, np1, iw)*twe2
123 continue
if (mlat == 0) goto 117
i = imid
do 124 k=1, nt
do 124 np1=mp1, ndo1, 2
br(mp1, np1, k) = br(mp1, np1, k)+wb(i, np1, iw)*we(i, 2*mp1-1, k)*wts(i)
bi(mp1, np1, k) = bi(mp1, np1, k)-wb(i, np1, iw)*we(i, 2*mp1-2, k)*wts(i)
124 continue
117 if (mp2 > ndo2) goto 120
do 121 k=1, nt
do 121 i=1, imm1
tve1 = ve(i, 2*mp1-1, k)*wts(i)
tve2 = ve(i, 2*mp1-2, k)*wts(i)
two1 = wo(i, 2*mp1-1, k)*wts(i)
two2 = wo(i, 2*mp1-2, k)*wts(i)
do 121 np1=mp2, ndo2, 2
br(mp1, np1, k) = br(mp1, np1, k)+vb(i, np1, iv)*tve2 &
                             +wb(i, np1, iw)*two1
bi(mp1, np1, k) = bi(mp1, np1, k)+vb(i, np1, iv)*tve1 &
                             -wb(i, np1, iw)*two2
121 continue
if (mlat == 0) goto 120
i = imid
do 122 k=1, nt
do 122 np1=mp2, ndo2, 2
br(mp1, np1, k) = br(mp1, np1, k)+vb(i, np1, iv)*ve(i, 2*mp1-2, k)*wts(i)
bi(mp1, np1, k) = bi(mp1, np1, k)+vb(i, np1, iv)*ve(i, 2*mp1-1, k)*wts(i)
122 continue
120 continue
return
!
!     case ityp=2 ,  no symmetries but br and bi equal zero   
!
200 call sphere_aux%vbin(0, nlat, nlon, 0, vb, iv, wvbin)
!
!     case m=0
!
do 215 k=1, nt
do 215 i=1, imid
tw = we(i, 1, k)*wts(i)
do 215 np1=2, ndo2, 2
cr(1, np1, k) = cr(1, np1, k)-vb(i, np1, iv)*tw
215 continue
do 216 k=1, nt
do 216 i=1, imm1
tw = wo(i, 1, k)*wts(i)
do 216 np1=3, ndo1, 2
cr(1, np1, k) = cr(1, np1, k)-vb(i, np1, iv)*tw
216 continue
!
!     case m = 1 through nlat-1
!
if (mmax < 2) return
do 220 mp1=2, mmax
m = mp1-1
mp2 = mp1+1
call sphere_aux%vbin(0, nlat, nlon, m, vb, iv, wvbin)
call sphere_aux%wbin(0, nlat, nlon, m, wb, iw, wwbin)
if (mp1 > ndo1) goto 217
do 223 k=1, nt
do 223 i=1, imm1
tve1 = ve(i, 2*mp1-1, k)*wts(i)
tve2 = ve(i, 2*mp1-2, k)*wts(i)
two1 = wo(i, 2*mp1-1, k)*wts(i)
two2 = wo(i, 2*mp1-2, k)*wts(i)
do 223 np1=mp1, ndo1, 2
cr(mp1, np1, k) = cr(mp1, np1, k)-vb(i, np1, iv)*two2 &
                             +wb(i, np1, iw)*tve1
ci(mp1, np1, k) = ci(mp1, np1, k)-vb(i, np1, iv)*two1 &
                             -wb(i, np1, iw)*tve2
223 continue
if (mlat == 0) goto 217
i = imid
do 224 k=1, nt
do 224 np1=mp1, ndo1, 2
cr(mp1, np1, k) = cr(mp1, np1, k)+wb(i, np1, iw)*ve(i, 2*mp1-1, k)*wts(i)
ci(mp1, np1, k) = ci(mp1, np1, k)-wb(i, np1, iw)*ve(i, 2*mp1-2, k)*wts(i)
224 continue
217 if (mp2 > ndo2) goto 220
do 221 k=1, nt
do 221 i=1, imm1
twe1 = we(i, 2*mp1-1, k)*wts(i)
twe2 = we(i, 2*mp1-2, k)*wts(i)
tvo1 = vo(i, 2*mp1-1, k)*wts(i)
tvo2 = vo(i, 2*mp1-2, k)*wts(i)
do 221 np1=mp2, ndo2, 2
cr(mp1, np1, k) = cr(mp1, np1, k)-vb(i, np1, iv)*twe2 &
                             +wb(i, np1, iw)*tvo1
ci(mp1, np1, k) = ci(mp1, np1, k)-vb(i, np1, iv)*twe1 &
                             -wb(i, np1, iw)*tvo2
221 continue
if (mlat == 0) goto 220
i = imid
do 222 k=1, nt
do 222 np1=mp2, ndo2, 2
cr(mp1, np1, k) = cr(mp1, np1, k)-vb(i, np1, iv)*we(i, 2*mp1-2, k)*wts(i)
ci(mp1, np1, k) = ci(mp1, np1, k)-vb(i, np1, iv)*we(i, 2*mp1-1, k)*wts(i)
222 continue
220 continue
return
!
!     case ityp=3 ,  v even , w odd
!
300 call sphere_aux%vbin(0, nlat, nlon, 0, vb, iv, wvbin)
!
!     case m=0
!
do 315 k=1, nt
do 315 i=1, imid
tv = ve(i, 1, k)*wts(i)
do 315 np1=2, ndo2, 2
br(1, np1, k) = br(1, np1, k)+vb(i, np1, iv)*tv
315 continue
do 316 k=1, nt
do 316 i=1, imm1
tw = wo(i, 1, k)*wts(i)
do 316 np1=3, ndo1, 2
cr(1, np1, k) = cr(1, np1, k)-vb(i, np1, iv)*tw
316 continue
!
!     case m = 1 through nlat-1
!
if (mmax < 2) return
do 320 mp1=2, mmax
m = mp1-1
mp2 = mp1+1
call sphere_aux%vbin(0, nlat, nlon, m, vb, iv, wvbin)
call sphere_aux%wbin(0, nlat, nlon, m, wb, iw, wwbin)
if (mp1 > ndo1) goto 317
do 323 k=1, nt
do 323 i=1, imm1
two1 = wo(i, 2*mp1-1, k)*wts(i)
two2 = wo(i, 2*mp1-2, k)*wts(i)
tve1 = ve(i, 2*mp1-1, k)*wts(i)
tve2 = ve(i, 2*mp1-2, k)*wts(i)
do 323 np1=mp1, ndo1, 2
cr(mp1, np1, k) = cr(mp1, np1, k)-vb(i, np1, iv)*two2 &
                             +wb(i, np1, iw)*tve1
ci(mp1, np1, k) = ci(mp1, np1, k)-vb(i, np1, iv)*two1 &
                             -wb(i, np1, iw)*tve2
323 continue
if (mlat == 0) goto 317
i = imid
do 324 k=1, nt
do 324 np1=mp1, ndo1, 2
cr(mp1, np1, k) = cr(mp1, np1, k)+wb(i, np1, iw)*ve(i, 2*mp1-1, k)*wts(i)
ci(mp1, np1, k) = ci(mp1, np1, k)-wb(i, np1, iw)*ve(i, 2*mp1-2, k)*wts(i)
324 continue
317 if (mp2 > ndo2) goto 320
do 321 k=1, nt
do 321 i=1, imm1
two1 = wo(i, 2*mp1-1, k)*wts(i)
two2 = wo(i, 2*mp1-2, k)*wts(i)
tve1 = ve(i, 2*mp1-1, k)*wts(i)
tve2 = ve(i, 2*mp1-2, k)*wts(i)
do 321 np1=mp2, ndo2, 2
br(mp1, np1, k) = br(mp1, np1, k)+vb(i, np1, iv)*tve2 &
                             +wb(i, np1, iw)*two1
bi(mp1, np1, k) = bi(mp1, np1, k)+vb(i, np1, iv)*tve1 &
                             -wb(i, np1, iw)*two2
321 continue
if (mlat == 0) goto 320
i = imid
do 322 k=1, nt
do 322 np1=mp2, ndo2, 2
br(mp1, np1, k) = br(mp1, np1, k)+vb(i, np1, iv)*ve(i, 2*mp1-2, k)*wts(i)
bi(mp1, np1, k) = bi(mp1, np1, k)+vb(i, np1, iv)*ve(i, 2*mp1-1, k)*wts(i)
322 continue
320 continue
return
!
!     case ityp=4 ,  v even, w odd, and cr and ci equal 0. 
!
400 call sphere_aux%vbin(1, nlat, nlon, 0, vb, iv, wvbin)
!
!     case m=0
!
do 415 k=1, nt
do 415 i=1, imid
tv = ve(i, 1, k)*wts(i)
do 415 np1=2, ndo2, 2
br(1, np1, k) = br(1, np1, k)+vb(i, np1, iv)*tv
415 continue
!
!     case m = 1 through nlat-1
!
if (mmax < 2) return
do 420 mp1=2, mmax
m = mp1-1
mp2 = mp1+1
call sphere_aux%vbin(1, nlat, nlon, m, vb, iv, wvbin)
call sphere_aux%wbin(1, nlat, nlon, m, wb, iw, wwbin)
if (mp2 > ndo2) goto 420
do 421 k=1, nt
do 421 i=1, imm1
two1 = wo(i, 2*mp1-1, k)*wts(i)
two2 = wo(i, 2*mp1-2, k)*wts(i)
tve1 = ve(i, 2*mp1-1, k)*wts(i)
tve2 = ve(i, 2*mp1-2, k)*wts(i)
do 421 np1=mp2, ndo2, 2
br(mp1, np1, k) = br(mp1, np1, k)+vb(i, np1, iv)*tve2 &
                             +wb(i, np1, iw)*two1
bi(mp1, np1, k) = bi(mp1, np1, k)+vb(i, np1, iv)*tve1 &
                             -wb(i, np1, iw)*two2
421 continue
if (mlat == 0) goto 420
i = imid
do 422 k=1, nt
do 422 np1=mp2, ndo2, 2
br(mp1, np1, k) = br(mp1, np1, k)+vb(i, np1, iv)*ve(i, 2*mp1-2, k)*wts(i)
bi(mp1, np1, k) = bi(mp1, np1, k)+vb(i, np1, iv)*ve(i, 2*mp1-1, k)*wts(i)
422 continue
420 continue
return
!
!     case ityp=5   v even, w odd, and br and bi equal zero
!
500 call sphere_aux%vbin(2, nlat, nlon, 0, vb, iv, wvbin)
!
!     case m=0
!
do 516 k=1, nt
do 516 i=1, imm1
tw = wo(i, 1, k)*wts(i)
do 516 np1=3, ndo1, 2
cr(1, np1, k) = cr(1, np1, k)-vb(i, np1, iv)*tw
516 continue
!
!     case m = 1 through nlat-1
!
if (mmax < 2) return
do 520 mp1=2, mmax
m = mp1-1
mp2 = mp1+1
call sphere_aux%vbin(2, nlat, nlon, m, vb, iv, wvbin)
call sphere_aux%wbin(2, nlat, nlon, m, wb, iw, wwbin)
if (mp1 > ndo1) goto 520
do 523 k=1, nt
do 523 i=1, imm1
two1 = wo(i, 2*mp1-1, k)*wts(i)
two2 = wo(i, 2*mp1-2, k)*wts(i)
tve1 = ve(i, 2*mp1-1, k)*wts(i)
tve2 = ve(i, 2*mp1-2, k)*wts(i)
do 523 np1=mp1, ndo1, 2
cr(mp1, np1, k) = cr(mp1, np1, k)-vb(i, np1, iv)*two2 &
                             +wb(i, np1, iw)*tve1
ci(mp1, np1, k) = ci(mp1, np1, k)-vb(i, np1, iv)*two1 &
                             -wb(i, np1, iw)*tve2
523 continue
if (mlat == 0) goto 520
i = imid
do 524 k=1, nt
do 524 np1=mp1, ndo1, 2
cr(mp1, np1, k) = cr(mp1, np1, k)+wb(i, np1, iw)*ve(i, 2*mp1-1, k)*wts(i)
ci(mp1, np1, k) = ci(mp1, np1, k)-wb(i, np1, iw)*ve(i, 2*mp1-2, k)*wts(i)
524 continue
520 continue
return
!
!     case ityp=6 ,  v odd , w even
!
600 call sphere_aux%vbin(0, nlat, nlon, 0, vb, iv, wvbin)
!
!     case m=0
!
do 615 k=1, nt
do 615 i=1, imid
tw = we(i, 1, k)*wts(i)
do 615 np1=2, ndo2, 2
cr(1, np1, k) = cr(1, np1, k)-vb(i, np1, iv)*tw
615 continue
do 616 k=1, nt
do 616 i=1, imm1
tv = vo(i, 1, k)*wts(i)
do 616 np1=3, ndo1, 2
br(1, np1, k) = br(1, np1, k)+vb(i, np1, iv)*tv
616 continue
!
!     case m = 1 through nlat-1
!
if (mmax < 2) return
do 620 mp1=2, mmax
m = mp1-1
mp2 = mp1+1
call sphere_aux%vbin(0, nlat, nlon, m, vb, iv, wvbin)
call sphere_aux%wbin(0, nlat, nlon, m, wb, iw, wwbin)
if (mp1 > ndo1) goto 617
do 623 k=1, nt
do 623 i=1, imm1
twe1 = we(i, 2*mp1-1, k)*wts(i)
twe2 = we(i, 2*mp1-2, k)*wts(i)
tvo1 = vo(i, 2*mp1-1, k)*wts(i)
tvo2 = vo(i, 2*mp1-2, k)*wts(i)
do 623 np1=mp1, ndo1, 2
br(mp1, np1, k) = br(mp1, np1, k)+vb(i, np1, iv)*tvo2 &
                             +wb(i, np1, iw)*twe1
bi(mp1, np1, k) = bi(mp1, np1, k)+vb(i, np1, iv)*tvo1 &
                             -wb(i, np1, iw)*twe2
623 continue
if (mlat == 0) goto 617
i = imid
do 624 k=1, nt
do 624 np1=mp1, ndo1, 2
br(mp1, np1, k) = br(mp1, np1, k)+wb(i, np1, iw)*we(i, 2*mp1-1, k)*wts(i)
bi(mp1, np1, k) = bi(mp1, np1, k)-wb(i, np1, iw)*we(i, 2*mp1-2, k)*wts(i)
624 continue
617 if (mp2 > ndo2) goto 620
do 621 k=1, nt
do 621 i=1, imm1
twe1 = we(i, 2*mp1-1, k)*wts(i)
twe2 = we(i, 2*mp1-2, k)*wts(i)
tvo1 = vo(i, 2*mp1-1, k)*wts(i)
tvo2 = vo(i, 2*mp1-2, k)*wts(i)
do 621 np1=mp2, ndo2, 2
cr(mp1, np1, k) = cr(mp1, np1, k)-vb(i, np1, iv)*twe2 &
                             +wb(i, np1, iw)*tvo1
ci(mp1, np1, k) = ci(mp1, np1, k)-vb(i, np1, iv)*twe1 &
                             -wb(i, np1, iw)*tvo2
621 continue
if (mlat == 0) goto 620
i = imid
do 622 k=1, nt
do 622 np1=mp2, ndo2, 2
cr(mp1, np1, k) = cr(mp1, np1, k)-vb(i, np1, iv)*we(i, 2*mp1-2, k)*wts(i)
ci(mp1, np1, k) = ci(mp1, np1, k)-vb(i, np1, iv)*we(i, 2*mp1-1, k)*wts(i)
622 continue
620 continue
return
!
!     case ityp=7   v odd, w even, and cr and ci equal zero
!
700 call sphere_aux%vbin(2, nlat, nlon, 0, vb, iv, wvbin)
!
!     case m=0
!
do 716 k=1, nt
do 716 i=1, imm1
tv = vo(i, 1, k)*wts(i)
do 716 np1=3, ndo1, 2
br(1, np1, k) = br(1, np1, k)+vb(i, np1, iv)*tv
716 continue
!
!     case m = 1 through nlat-1
!
if (mmax < 2) return
do 720 mp1=2, mmax
m = mp1-1
mp2 = mp1+1
call sphere_aux%vbin(2, nlat, nlon, m, vb, iv, wvbin)
call sphere_aux%wbin(2, nlat, nlon, m, wb, iw, wwbin)
if (mp1 > ndo1) goto 720
do 723 k=1, nt
do 723 i=1, imm1
twe1 = we(i, 2*mp1-1, k)*wts(i)
twe2 = we(i, 2*mp1-2, k)*wts(i)
tvo1 = vo(i, 2*mp1-1, k)*wts(i)
tvo2 = vo(i, 2*mp1-2, k)*wts(i)
do 723 np1=mp1, ndo1, 2
br(mp1, np1, k) = br(mp1, np1, k)+vb(i, np1, iv)*tvo2 &
                             +wb(i, np1, iw)*twe1
bi(mp1, np1, k) = bi(mp1, np1, k)+vb(i, np1, iv)*tvo1 &
                             -wb(i, np1, iw)*twe2
723 continue
if (mlat == 0) goto 720
i = imid
do 724 k=1, nt
do 724 np1=mp1, ndo1, 2
br(mp1, np1, k) = br(mp1, np1, k)+wb(i, np1, iw)*we(i, 2*mp1-1, k)*wts(i)
bi(mp1, np1, k) = bi(mp1, np1, k)-wb(i, np1, iw)*we(i, 2*mp1-2, k)*wts(i)
724 continue
720 continue
return
!
!     case ityp=8   v odd, w even, and both br and bi equal zero
!
800 call sphere_aux%vbin(1, nlat, nlon, 0, vb, iv, wvbin)
!
!     case m=0
!
do 815 k=1, nt
do 815 i=1, imid
tw = we(i, 1, k)*wts(i)
do 815 np1=2, ndo2, 2
cr(1, np1, k) = cr(1, np1, k)-vb(i, np1, iv)*tw
815 continue
!
!     case m = 1 through nlat-1
!
if (mmax < 2) return
do 820 mp1=2, mmax
m = mp1-1
mp2 = mp1+1
call sphere_aux%vbin(1, nlat, nlon, m, vb, iv, wvbin)
call sphere_aux%wbin(1, nlat, nlon, m, wb, iw, wwbin)
if (mp2 > ndo2) goto 820
do 821 k=1, nt
do 821 i=1, imm1
twe1 = we(i, 2*mp1-1, k)*wts(i)
twe2 = we(i, 2*mp1-2, k)*wts(i)
tvo1 = vo(i, 2*mp1-1, k)*wts(i)
tvo2 = vo(i, 2*mp1-2, k)*wts(i)
do 821 np1=mp2, ndo2, 2
cr(mp1, np1, k) = cr(mp1, np1, k)-vb(i, np1, iv)*twe2 &
                             +wb(i, np1, iw)*tvo1
ci(mp1, np1, k) = ci(mp1, np1, k)-vb(i, np1, iv)*twe1 &
                             -wb(i, np1, iw)*tvo2
821 continue
if (mlat == 0) goto 820
i = imid
do 822 k=1, nt
do 822 np1=mp2, ndo2, 2
cr(mp1, np1, k) = cr(mp1, np1, k)-vb(i, np1, iv)*we(i, 2*mp1-2, k)*wts(i)
ci(mp1, np1, k) = ci(mp1, np1, k)-vb(i, np1, iv)*we(i, 2*mp1-1, k)*wts(i)
822 continue
820 continue

end subroutine vhagc1

end subroutine vhagc

subroutine vhagci(nlat, nlon, wvhagc, lvhagc, dwork, ldwork, ierror)

integer (ip) :: ierror
integer (ip) :: imid
integer (ip) :: iw1
integer (ip) :: iw2
integer (ip) :: iw3
integer (ip) :: iwrk
integer (ip) :: jw1
integer (ip) :: jw2
integer (ip) :: jw3
integer (ip) :: labc
integer (ip) :: ldwork
integer (ip) :: lvhagc
integer (ip) :: lwk
integer (ip) :: lwvbin
integer (ip) :: lzz1
integer (ip) :: mmax
integer (ip) :: nlat
integer (ip) :: nlon
real (wp) :: wvhagc(lvhagc)
real (wp) :: dwork(ldwork)
real (wp) :: dummy_variable

type (HFFTpack)      :: hfft
type (SpherepackAux) :: sphere_aux

ierror = 1
if (nlat < 3) return
ierror = 2
if (nlon < 1) return
ierror = 3
imid = (nlat+1)/2
lzz1 = 2*nlat*imid
mmax = min(nlat, (nlon+1)/2)
labc = 3*(max(mmax-2, 0)*(nlat+nlat-mmax-1))/2
imid = (nlat+1)/2
if (lvhagc < 2*(lzz1+labc)+nlon+imid+15) return
ierror = 4
if (ldwork < 2*nlat*(nlat+1)+1) return
ierror = 0
!
!     compute gaussian points in first nlat+1 words of dwork
!     real
!
lwk = nlat*(nlat+2)

jw1 = 1
!     jw2 = jw1+nlat+nlat
!     jw3 = jw2+nlat+nlat
jw2 = jw1+nlat
jw3 = jw2+nlat
call compute_gaussian_latitudes_and_weights(nlat, dwork(jw1), dwork(jw2), dummy_variable, lwk, ierror)
imid = (nlat+1)/2
!
!     set first imid words of real weights in dwork
!     as single precision in first imid words of wvhagc
!
call setwts(imid, dwork(nlat+1), wvhagc)
!
!     first nlat+1 words of dwork contain  double theta
!
!     iwrk = nlat+2
iwrk = (nlat+1)/2 +1
iw1 = imid+1
lwvbin = lzz1+labc
iw2 = iw1+lwvbin
iw3 = iw2+lwvbin

call sphere_aux%vbgint(nlat, nlon, dwork, wvhagc(iw1), dwork(iwrk))

call sphere_aux%wbgint(nlat, nlon, dwork, wvhagc(iw2), dwork(iwrk))

call hfft%initialize(nlon, wvhagc(iw3))

end subroutine vhagci



subroutine setwts(imid, dwts, wts)

integer (ip), intent (in)  :: imid
real (wp),    intent (in)  :: dwts(imid)
real (wp),    intent (out) :: wts(imid)
!
!     set first imid =(nlat+1)/2 of real weights in dwts
!     as single precision in wts
!

wts(1:imid) = dwts(1:imid)

end subroutine setwts

end module module_vhagc
