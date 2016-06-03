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
! ... file vhags.f
!
!     this file contains code and documentation for subroutines
!     vhags and vhagsi
!
! ... files which must be loaded with vhags.f
!
!     type_SpherepackAux.f, type_HFFTpack.f, gaqd.f
!
!     subroutine vhags(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci,
!    +                 mdab, ndab, wvhags, lvhags, work, lwork, ierror)
!
!     subroutine vhags performs the vector spherical harmonic analysis
!     on the vector field (v, w) and stores the result in the arrays
!     br, bi, cr, and ci. v(i, j) and w(i, j) are the colatitudinal
!     (measured from the north pole) and east longitudinal components
!     respectively, located at the gaussian colatitude point theta(i)
!     and longitude phi(j) = (j-1)*2*pi/nlon. the spectral
!     representation of (v, w) is given at output parameters v, w in
!     subroutine vhses.  
!
!     input parameters
!
!     nlat   the number of points in the gaussian colatitude grid on the
!            full sphere. these lie in the interval (0, pi) and are computed
!            in radians in theta(1) <...< theta(nlat) by subroutine gaqd.
!            if nlat is odd the equator will be included as the grid point
!            theta((nlat+1)/2).  if nlat is even the equator will be
!            excluded as a grid point and will lie half way between
!            theta(nlat/2) and theta(nlat/2+1). nlat must be at least 3.
!            note: on the half sphere, the number of grid points in the
!            colatitudinal direction is nlat/2 if nlat is even or
!            (nlat+1)/2 if nlat is odd.
!
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
!     nt     the number of analyses.  in the program that calls vhags,
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
!            components at the gaussian colatitude point theta(i)
!            and longitude phi(j) = (j-1)*2*pi/nlon. the index ranges
!            are defined above at the input parameter ityp.
!
!     idvw   the first dimension of the arrays v, w as it appears in
!            the program that calls vhags. if ityp .le. 2 then idvw
!            must be at least nlat.  if ityp .gt. 2 and nlat is
!            even then idvw must be at least nlat/2. if ityp .gt. 2
!            and nlat is odd then idvw must be at least (nlat+1)/2.
!
!     jdvw   the second dimension of the arrays v, w as it appears in
!            the program that calls vhags. jdvw must be at least nlon.
!
!     mdab   the first dimension of the arrays br, bi, cr, and ci as it
!            appears in the program that calls vhags. mdab must be at
!            least min(nlat, nlon/2) if nlon is even or at least
!            min(nlat, (nlon+1)/2) if nlon is odd.
!
!     ndab   the second dimension of the arrays br, bi, cr, and ci as it
!            appears in the program that calls vhags. ndab must be at
!            least nlat.
!
!     wvhags an array which must be initialized by subroutine vhgsi.
!            once initialized, wvhags can be used repeatedly by vhags
!            as long as nlon and nlat remain unchanged.  wvhags must
!            not be altered between calls of vhags.
!
!     lvhags the dimension of the array wvhags as it appears in the
!            program that calls vhags. define
!
!               l1 = min(nlat, nlon/2) if nlon is even or
!               l1 = min(nlat, (nlon+1)/2) if nlon is odd
!
!            and
!
!               l2 = nlat/2        if nlat is even or
!               l2 = (nlat+1)/2    if nlat is odd
!
!            then lvhags must be at least
!
!            l1*l2(nlat+nlat-l1+1)+nlon+15
!
!        ??? (nlat+1)*(nlat+1)*nlat/2 + nlon + 15
!
!
!
!     work   a work array that does not have to be saved.
!
!     lwork  the dimension of the array work as it appears in the
!            program that calls vhags. define
!
!               l2 = nlat/2        if nlat is even or
!               l2 = (nlat+1)/2    if nlat is odd
!
!            if ityp .le. 2 then lwork must be at least
!            the larger of the two quantities
!
!               3*nlat*(nlat+1)+2  (required by vhagsi)
!
!            and
!
!               (2*nt+1)*nlat*nlon
!
!            if ityp .gt. 2 then lwork must be at least
!            the larger of the two quantities
!
!               3*nlat*(nlat+1)+2  (required by vhagsi)
!
!            and
!
!              (2*nt+1)*l2*nlon
!
!
!     **************************************************************
!
!     output parameters
!
!     br, bi  two or three dimensional arrays (see input parameter nt)
!     cr, ci  that contain the vector spherical harmonic coefficients
!            in the spectral representation of v(i, j) and w(i, j) given
!            in the discription of subroutine vhses. br(mp1, np1),
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
!            = 9  error in the specification of lvhags
!            = 10 error in the specification of lwork
!
!
!     subroutine vhagsi(nlat, nlon, wvhags, lvhags, work, lwork, ierror)
!
!     subroutine vhagsi initializes the array wvhags which can then be
!     used repeatedly by subroutine vhags until nlat or nlon is changed.
!
!     input parameters
!
!     nlat   the number of points in the gaussian colatitude grid on the
!            full sphere. these lie in the interval (0, pi) and are computed
!            in radians in theta(1) <...< theta(nlat) by subroutine gaqd.
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
!     lvhags the dimension of the array wvhags as it appears in the
!            program that calls vhagsi.  define
!
!               l1 = min(nlat, nlon/2) if nlon is even or
!               l1 = min(nlat, (nlon+1)/2) if nlon is odd
!
!            and
!
!               l2 = nlat/2        if nlat is even or
!               l2 = (nlat+1)/2    if nlat is odd
!
!            then lvhags must be at least
!
!               3*nlat*(nlat+1)+2  (required by vhagsi)
!
!     dwork  a real work space that does not need to be saved
!
!     ldwork the dimension of the array dwork as it appears in the
!            program that calls vhagsi. ldwork must be at least
!
!                   (3*nlat*(nlat+3)+2)/2
!
!     **************************************************************
!
!     output parameters
!
!     wvhags an array which is initialized for use by subroutine vhags.
!            once initialized, wvhags can be used repeatedly by vhags
!            as long as nlat and nlon remain unchanged.  wvhags must not
!            be altered between calls of vhags.
!
!
!     ierror = 0  no errors
!            = 1  error in the specification of nlat
!            = 2  error in the specification of nlon
!            = 3  error in the specification of lvhags
!            = 4  error in the specification of ldwork
!
module module_vhags

    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32

    use type_HFFTpack, only: &
        HFFTpack

    use type_SpherepackAux, only: &
        SpherepackAux

    use module_gaqd, only: &
        gaqd

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    public :: vhags
    public :: vhagsi

contains

    subroutine vhags(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, &
        mdab, ndab, wvhags, lvhags, work, lwork, ierror)

        real (wp) :: bi
        real (wp) :: br
        real (wp) :: ci
        real (wp) :: cr
        integer (ip) :: idv
        integer (ip) :: idvw
        integer (ip) :: idz
        integer (ip) :: ierror
        integer (ip) :: imid
        integer (ip) :: ist
        integer (ip) :: ityp
        integer (ip) :: iw1
        integer (ip) :: iw2
        integer (ip) :: iw3
        integer (ip) :: iw4
        integer (ip) :: jdvw
        integer (ip) :: jw1
        integer (ip) :: jw2
        integer (ip) :: jw3
        integer (ip) :: lmn
        integer (ip) :: lnl
        integer (ip) :: lvhags
        integer (ip) :: lwork
        integer (ip) :: lzimn
        integer (ip) :: mdab
        integer (ip) :: mmax
        integer (ip) :: ndab
        integer (ip) :: nlat
        integer (ip) :: nlon
        integer (ip) :: nt
        real (wp) :: v
        real (wp) :: w
        real (wp) :: work
        real (wp) :: wvhags
        dimension v(idvw, jdvw, *), w(idvw, jdvw,*), br(mdab,ndab,*), &
            bi(mdab, ndab,*), cr(mdab, ndab,*), ci(mdab, ndab,*), &
            work(*), wvhags(*)
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
        idz = (mmax*(nlat+nlat-mmax+1))/2
        lzimn = idz*imid
        if (lvhags < lzimn+lzimn+nlon+15) return
        ierror = 10
        idv = nlat
        if (ityp > 2) idv = imid
        lnl = nt*idv*nlon
        if (lwork < lnl+lnl+idv*nlon) return
        ierror = 0

        if (ityp <= 2) then
            ist = imid
        else
            ist = 0
        end if

        !
        !     set wvhags pointers
        !
        lmn = nlat*(nlat+1)/2
        jw1 = 1
        jw2 = jw1+imid*lmn
        jw3 = jw2+imid*lmn
        !
        !     set work pointers
        !
        iw1 = ist+1
        iw2 = lnl+1
        iw3 = iw2+ist
        iw4 = iw2+lnl
        call vhags1(nlat, nlon, ityp, nt, imid, idvw, jdvw, v, w, mdab, ndab, &
            br, bi, cr, ci, idv, work, work(iw1), work(iw2), work(iw3), &
            work(iw4), idz, wvhags(jw1), wvhags(jw2), wvhags(jw3))

    contains

        subroutine vhags1(nlat, nlon, ityp, nt, imid, idvw, jdvw, v, w, mdab, &
            ndab, br, bi, cr, ci, idv, ve, vo, we, wo, work, idz, vb, wb, wrfft)

            real (wp) :: bi
            real (wp) :: br
            real (wp) :: ci
            real (wp) :: cr
            real (wp) :: fsn
            integer (ip) :: i
            integer (ip) :: idv
            integer (ip) :: idvw
            integer (ip) :: idz
            integer (ip) :: imid
            integer (ip) :: imm1
            integer (ip) :: ityp
            integer (ip) :: j
            integer (ip) :: jdvw
            integer (ip) :: k
            integer (ip) :: m
            integer (ip) :: mb
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
            real (wp) :: v
            real (wp) :: vb
            real (wp) :: ve
            real (wp) :: vo
            real (wp) :: w
            real (wp) :: wb
            real (wp) :: we
            real (wp) :: wo
            real (wp) :: work
            real (wp) :: wrfft
            dimension v(idvw, jdvw, *), w(idvw, jdvw, *), br(mdab, ndab,*), &
                bi(mdab, ndab, *), cr(mdab, ndab, *), ci(mdab, ndab, *), &
                ve(idv, nlon, *), vo(idv, nlon, *), we(idv, nlon, *), &
                wo(idv, nlon, *), work(*), &
                vb(imid,*), wb(imid,*), wrfft(*)

            type (HFFTpack)      :: hfft

            nlp1 = nlat+1
            tsn = 2.0/nlon
            fsn = 4.0/nlon
            mlat = mod(nlat, 2)
            mlon = mod(nlon, 2)
            mmax = min(nlat, (nlon+1)/2)

            select case (mlat)
                case (0)
                    imm1 = imid
                    ndo1 = nlat
                    ndo2 = nlat-1
                case default
                    imm1 = imid-1
                    ndo1 = nlat-1
                    ndo2 = nlat
            end select

            if (ityp <= 2) then
                do k=1, nt
                    do i=1, imm1
                        do j=1, nlon
                            ve(i, j, k) = tsn*(v(i, j, k)+v(nlp1-i, j, k))
                            vo(i, j, k) = tsn*(v(i, j, k)-v(nlp1-i, j, k))
                            we(i, j, k) = tsn*(w(i, j, k)+w(nlp1-i, j, k))
                            wo(i, j, k) = tsn*(w(i, j, k)-w(nlp1-i, j, k))
                        end do
                    end do
                end do
            else
                do k=1, nt
                    do i=1, imm1
                        do j=1, nlon
                            ve(i, j, k) = fsn*v(i, j, k)
                            vo(i, j, k) = fsn*v(i, j, k)
                            we(i, j, k) = fsn*w(i, j, k)
                            wo(i, j, k) = fsn*w(i, j, k)
                        end do
                    end do
                end do
            end if

            if (mlat /= 0) then
                do k=1, nt
                    do j=1, nlon
                        ve(imid, j, k) = tsn*v(imid, j, k)
                        we(imid, j, k) = tsn*w(imid, j, k)
                    end do
                end do
            end if

            do k=1, nt
                call hfft%forward(idv, nlon, ve(1, 1, k), idv, wrfft, work)
                call hfft%forward(idv, nlon, we(1, 1, k), idv, wrfft, work)
            end do

            !
            !==> Set polar coefficients to zero
            !
            select case (ityp)
                case (0:1,3:4,6:7)
                    do k=1, nt
                        do mp1=1, mmax
                            do np1=mp1, nlat
                                br(mp1, np1, k) = 0.0
                                bi(mp1, np1, k) = 0.0
                            end do
                        end do
                    end do
            end select

            !
            !==> Set azimuthal coefficients to zero
            !
            select case (ityp)
                case (0,2:3,5:6,8)
                    do k=1, nt
                        do mp1=1, mmax
                            do np1=mp1, nlat
                                cr(mp1, np1, k) = 0.0
                                ci(mp1, np1, k) = 0.0
                            end do
                        end do
                    end do
            end select

            select case (ityp)
                case (0)
                    !
                    !==> case ityp=0 ,  no symmetries
                    !
                    !    case m=0
                    !
                    do k=1, nt
                        do i=1, imid
                            do np1=2, ndo2, 2
                                br(1, np1, k) = br(1, np1, k)+vb(i, np1)*ve(i, 1, k)
                                cr(1, np1, k) = cr(1, np1, k)-vb(i, np1)*we(i, 1, k)
                            end do
                        end do
                    end do

                    do k=1, nt
                        do i=1, imm1
                            do np1=3, ndo1, 2
                                br(1, np1, k) = br(1, np1, k)+vb(i, np1)*vo(i, 1, k)
                                cr(1, np1, k) = cr(1, np1, k)-vb(i, np1)*wo(i, 1, k)
                            end do
                        end do
                    end do
                    !
                    !==> case m = 1 through nlat-1
                    !
                    if (mmax < 2) then
                        return
                    end if

                    do mp1=2, mmax
                        m = mp1-1
                        mb = m*nlat-(m*(m+1))/2
                        mp2 = mp1+1

                        if (mp1 <= ndo1) then

                            do k=1, nt
                                do i=1, imm1
                                    do np1=mp1, ndo1, 2
                                        br(mp1, np1, k) = br(mp1, np1, k)+vb(i, np1+mb)*vo(i, 2*mp1-2, k) &
                                            +wb(i, np1+mb)*we(i, 2*mp1-1, k)
                                        bi(mp1, np1, k) = bi(mp1, np1, k)+vb(i, np1+mb)*vo(i, 2*mp1-1, k) &
                                            -wb(i, np1+mb)*we(i, 2*mp1-2, k)
                                        cr(mp1, np1, k) = cr(mp1, np1, k)-vb(i, np1+mb)*wo(i, 2*mp1-2, k) &
                                            +wb(i, np1+mb)*ve(i, 2*mp1-1, k)
                                        ci(mp1, np1, k) = ci(mp1, np1, k)-vb(i, np1+mb)*wo(i, 2*mp1-1, k) &
                                            -wb(i, np1+mb)*ve(i, 2*mp1-2, k)
                                    end do
                                end do
                            end do

                            if (mlat /= 0) then
                                do k=1, nt
                                    do np1=mp1, ndo1, 2
                                        br(mp1, np1, k) = br(mp1, np1, k)+wb(imid, np1+mb)*we(imid, 2*mp1-1, k)
                                        bi(mp1, np1, k) = bi(mp1, np1, k)-wb(imid, np1+mb)*we(imid, 2*mp1-2, k)
                                        cr(mp1, np1, k) = cr(mp1, np1, k)+wb(imid, np1+mb)*ve(imid, 2*mp1-1, k)
                                        ci(mp1, np1, k) = ci(mp1, np1, k)-wb(imid, np1+mb)*ve(imid, 2*mp1-2, k)
                                    end do
                                end do

                            end if
                        end if

                        if (mp2 > ndo2) exit

                        do k=1, nt
                            do i=1, imm1
                                do np1=mp2, ndo2, 2
                                    br(mp1, np1, k) = br(mp1, np1, k)+vb(i, np1+mb)*ve(i, 2*mp1-2, k) &
                                        +wb(i, np1+mb)*wo(i, 2*mp1-1, k)
                                    bi(mp1, np1, k) = bi(mp1, np1, k)+vb(i, np1+mb)*ve(i, 2*mp1-1, k) &
                                        -wb(i, np1+mb)*wo(i, 2*mp1-2, k)
                                    cr(mp1, np1, k) = cr(mp1, np1, k)-vb(i, np1+mb)*we(i, 2*mp1-2, k) &
                                        +wb(i, np1+mb)*vo(i, 2*mp1-1, k)
                                    ci(mp1, np1, k) = ci(mp1, np1, k)-vb(i, np1+mb)*we(i, 2*mp1-1, k) &
                                        -wb(i, np1+mb)*vo(i, 2*mp1-2, k)
                                end do
                            end do
                        end do

                        if (mlat == 0) exit

                        do k=1, nt
                            do np1=mp2, ndo2, 2
                                br(mp1, np1, k) = br(mp1, np1, k)+vb(imid, np1+mb)*ve(imid, 2*mp1-2, k)
                                bi(mp1, np1, k) = bi(mp1, np1, k)+vb(imid, np1+mb)*ve(imid, 2*mp1-1, k)
                                cr(mp1, np1, k) = cr(mp1, np1, k)-vb(imid, np1+mb)*we(imid, 2*mp1-2, k)
                                ci(mp1, np1, k) = ci(mp1, np1, k)-vb(imid, np1+mb)*we(imid, 2*mp1-1, k)
                            end do
                        end do
                    end do
                case (1)
                    !
                    !==> case ityp=1 ,  no symmetries but cr and ci equal zero
                    !
                    !    case m=0
                    !
                    do k=1, nt
                        do i=1, imid
                            do np1=2, ndo2, 2
                                br(1, np1, k) = br(1, np1, k)+vb(i, np1)*ve(i, 1, k)
                            end do
                        end do
                    end do

                    do k=1, nt
                        do i=1, imm1
                            do np1=3, ndo1, 2
                                br(1, np1, k) = br(1, np1, k)+vb(i, np1)*vo(i, 1, k)
                            end do
                        end do
                    end do
                    !
                    !==> case m = 1 through nlat-1
                    !
                    if (mmax < 2) return

                    do mp1=2, mmax
                        m = mp1-1
                        mb = m*nlat-(m*(m+1))/2
                        mp2 = mp1+1

                        if (mp1 <= ndo1) then

                            do k=1, nt
                                do i=1, imm1
                                    do np1=mp1, ndo1, 2
                                        br(mp1, np1, k) = br(mp1, np1, k)+vb(i, np1+mb)*vo(i, 2*mp1-2, k) &
                                            +wb(i, np1+mb)*we(i, 2*mp1-1, k)
                                        bi(mp1, np1, k) = bi(mp1, np1, k)+vb(i, np1+mb)*vo(i, 2*mp1-1, k) &
                                            -wb(i, np1+mb)*we(i, 2*mp1-2, k)
                                    end do
                                end do
                            end do

                            if (mlat /= 0) then
                                do k=1, nt
                                    do np1=mp1, ndo1, 2
                                        br(mp1, np1, k) = br(mp1, np1, k)+wb(imid, np1+mb)*we(imid, 2*mp1-1, k)
                                        bi(mp1, np1, k) = bi(mp1, np1, k)-wb(imid, np1+mb)*we(imid, 2*mp1-2, k)
                                    end do
                                end do
                            end if
                        end if

                        if (mp2 > ndo2) exit

                        do k=1, nt
                            do i=1, imm1
                                do np1=mp2, ndo2, 2
                                    br(mp1, np1, k) = br(mp1, np1, k)+vb(i, np1+mb)*ve(i, 2*mp1-2, k) &
                                        +wb(i, np1+mb)*wo(i, 2*mp1-1, k)
                                    bi(mp1, np1, k) = bi(mp1, np1, k)+vb(i, np1+mb)*ve(i, 2*mp1-1, k) &
                                        -wb(i, np1+mb)*wo(i, 2*mp1-2, k)
                                end do
                            end do
                        end do

                        if (mlat == 0) exit

                        do k=1, nt
                            do np1=mp2, ndo2, 2
                                br(mp1, np1, k) = br(mp1, np1, k)+vb(imid, np1+mb)*ve(imid, 2*mp1-2, k)
                                bi(mp1, np1, k) = bi(mp1, np1, k)+vb(imid, np1+mb)*ve(imid, 2*mp1-1, k)
                            end do
                        end do
                    end do
                case (2)
                    !
                    !==> case ityp=2 ,  no symmetries but br and bi equal zero
                    !
                    !    case m=0
                    !
                    do k=1, nt
                        do i=1, imid
                            do np1=2, ndo2, 2
                                cr(1, np1, k) = cr(1, np1, k)-vb(i, np1)*we(i, 1, k)
                            end do
                        end do
                    end do

                    do k=1, nt
                        do i=1, imm1
                            do np1=3, ndo1, 2
                                cr(1, np1, k) = cr(1, np1, k)-vb(i, np1)*wo(i, 1, k)
                            end do
                        end do
                    end do
                    !
                    !==> case m = 1 through nlat-1
                    !
                    if (mmax < 2) return

                    do mp1=2, mmax
                        m = mp1-1
                        mb = m*nlat-(m*(m+1))/2
                        mp2 = mp1+1

                        if (mp1 <= ndo1) then

                            do k=1, nt
                                do i=1, imm1
                                    do np1=mp1, ndo1, 2
                                        cr(mp1, np1, k) = cr(mp1, np1, k)-vb(i, np1+mb)*wo(i, 2*mp1-2, k) &
                                            +wb(i, np1+mb)*ve(i, 2*mp1-1, k)
                                        ci(mp1, np1, k) = ci(mp1, np1, k)-vb(i, np1+mb)*wo(i, 2*mp1-1, k) &
                                            -wb(i, np1+mb)*ve(i, 2*mp1-2, k)
                                    end do
                                end do
                            end do

                            if (mlat /= 0) then
                                do k=1, nt
                                    do np1=mp1, ndo1, 2
                                        cr(mp1, np1, k) = cr(mp1, np1, k)+wb(imid, np1+mb)*ve(imid, 2*mp1-1, k)
                                        ci(mp1, np1, k) = ci(mp1, np1, k)-wb(imid, np1+mb)*ve(imid, 2*mp1-2, k)
                                    end do
                                end do
                            end if
                        end if

                        if (mp2 > ndo2) exit

                        do k=1, nt
                            do i=1, imm1
                                do np1=mp2, ndo2, 2
                                    cr(mp1, np1, k) = cr(mp1, np1, k)-vb(i, np1+mb)*we(i, 2*mp1-2, k) &
                                        +wb(i, np1+mb)*vo(i, 2*mp1-1, k)
                                    ci(mp1, np1, k) = ci(mp1, np1, k)-vb(i, np1+mb)*we(i, 2*mp1-1, k) &
                                        -wb(i, np1+mb)*vo(i, 2*mp1-2, k)
                                end do
                            end do
                        end do

                        if (mlat == 0) exit

                        do k=1, nt
                            do np1=mp2, ndo2, 2
                                cr(mp1, np1, k) = cr(mp1, np1, k)-vb(imid, np1+mb)*we(imid, 2*mp1-2, k)
                                ci(mp1, np1, k) = ci(mp1, np1, k)-vb(imid, np1+mb)*we(imid, 2*mp1-1, k)
                            end do
                        end do
                    end do
                case (3)
                    !
                    !==> case ityp=3 ,  v even , w odd
                    !
                    !    case m=0
                    !
                    do k=1, nt
                        do i=1, imid
                            do np1=2, ndo2, 2
                                br(1, np1, k) = br(1, np1, k)+vb(i, np1)*ve(i, 1, k)
                            end do
                        end do
                    end do

                    do k=1, nt
                        do i=1, imm1
                            do np1=3, ndo1, 2
                                cr(1, np1, k) = cr(1, np1, k)-vb(i, np1)*wo(i, 1, k)
                            end do
                        end do
                    end do
                    !
                    !==> case m = 1 through nlat-1
                    !
                    if (mmax < 2) return

                    do mp1=2, mmax
                        m = mp1-1
                        mb = m*nlat-(m*(m+1))/2
                        mp2 = mp1+1

                        if (mp1 <= ndo1) then

                            do k=1, nt
                                do i=1, imm1
                                    do np1=mp1, ndo1, 2
                                        cr(mp1, np1, k) = cr(mp1, np1, k)-vb(i, np1+mb)*wo(i, 2*mp1-2, k) &
                                            +wb(i, np1+mb)*ve(i, 2*mp1-1, k)
                                        ci(mp1, np1, k) = ci(mp1, np1, k)-vb(i, np1+mb)*wo(i, 2*mp1-1, k) &
                                            -wb(i, np1+mb)*ve(i, 2*mp1-2, k)
                                    end do
                                end do
                            end do

                            if (mlat /= 0) then

                                do k=1, nt
                                    do np1=mp1, ndo1, 2
                                        cr(mp1, np1, k) = cr(mp1, np1, k)+wb(imid, np1+mb)*ve(imid, 2*mp1-1, k)
                                        ci(mp1, np1, k) = ci(mp1, np1, k)-wb(imid, np1+mb)*ve(imid, 2*mp1-2, k)
                                    end do
                                end do
                            end if
                        end if

                        if (mp2 > ndo2) exit

                        do k=1, nt
                            do i=1, imm1
                                do np1=mp2, ndo2, 2
                                    br(mp1, np1, k) = br(mp1, np1, k)+vb(i, np1+mb)*ve(i, 2*mp1-2, k) &
                                        +wb(i, np1+mb)*wo(i, 2*mp1-1, k)
                                    bi(mp1, np1, k) = bi(mp1, np1, k)+vb(i, np1+mb)*ve(i, 2*mp1-1, k) &
                                        -wb(i, np1+mb)*wo(i, 2*mp1-2, k)
                                end do
                            end do
                        end do

                        if (mlat == 0) exit

                        do k=1, nt
                            do np1=mp2, ndo2, 2
                                br(mp1, np1, k) = br(mp1, np1, k)+vb(imid, np1+mb)*ve(imid, 2*mp1-2, k)
                                bi(mp1, np1, k) = bi(mp1, np1, k)+vb(imid, np1+mb)*ve(imid, 2*mp1-1, k)
                            end do
                        end do
                    end do
                case (4)
                    !
                    !==> case ityp=4 ,  v even, w odd, and cr and ci equal 0.
                    !
                    !     case m=0
                    !
                    do k=1, nt
                        do i=1, imid
                            do np1=2, ndo2, 2
                                br(1, np1, k) = br(1, np1, k)+vb(i, np1)*ve(i, 1, k)
                            end do
                        end do
                    end do
                    !
                    !==> case m = 1 through nlat-1
                    !
                    if (mmax < 2) return

                    do mp1=2, mmax
                        m = mp1-1
                        mb = m*nlat-(m*(m+1))/2
                        mp2 = mp1+1

                        if (mp2 > ndo2) exit

                        do k=1, nt
                            do i=1, imm1
                                do np1=mp2, ndo2, 2
                                    br(mp1, np1, k) = br(mp1, np1, k)+vb(i, np1+mb)*ve(i, 2*mp1-2, k) &
                                        +wb(i, np1+mb)*wo(i, 2*mp1-1, k)
                                    bi(mp1, np1, k) = bi(mp1, np1, k)+vb(i, np1+mb)*ve(i, 2*mp1-1, k) &
                                        -wb(i, np1+mb)*wo(i, 2*mp1-2, k)
                                end do
                            end do
                        end do

                        if (mlat == 0) exit

                        do k=1, nt
                            do np1=mp2, ndo2, 2
                                br(mp1, np1, k) = br(mp1, np1, k)+vb(imid, np1+mb)*ve(imid, 2*mp1-2, k)
                                bi(mp1, np1, k) = bi(mp1, np1, k)+vb(imid, np1+mb)*ve(imid, 2*mp1-1, k)
                            end do
                        end do
                    end do
                case (5)
                    !
                    !==> case ityp=5   v even, w odd, and br and bi equal zero
                    !
                    !     case m=0
                    !
                    do k=1, nt
                        do i=1, imm1
                            do np1=3, ndo1, 2
                                cr(1, np1, k) = cr(1, np1, k)-vb(i, np1)*wo(i, 1, k)
                            end do
                        end do
                    end do
                    !
                    !==> case m = 1 through nlat-1
                    !
                    if (mmax < 2) return

                    do mp1=2, mmax
                        m = mp1-1
                        mb = m*nlat-(m*(m+1))/2
                        mp2 = mp1+1

                        if (mp1 > ndo1) exit

                        do  k=1, nt
                            do  i=1, imm1
                                do  np1=mp1, ndo1, 2
                                    cr(mp1, np1, k) = cr(mp1, np1, k)-vb(i, np1+mb)*wo(i, 2*mp1-2, k) &
                                        +wb(i, np1+mb)*ve(i, 2*mp1-1, k)
                                    ci(mp1, np1, k) = ci(mp1, np1, k)-vb(i, np1+mb)*wo(i, 2*mp1-1, k) &
                                        -wb(i, np1+mb)*ve(i, 2*mp1-2, k)
                                end do
                            end do
                        end do

                        if (mlat == 0) exit

                        do k=1, nt
                            do np1=mp1, ndo1, 2
                                cr(mp1, np1, k) = cr(mp1, np1, k)+wb(imid, np1+mb)*ve(imid, 2*mp1-1, k)
                                ci(mp1, np1, k) = ci(mp1, np1, k)-wb(imid, np1+mb)*ve(imid, 2*mp1-2, k)
                            end do
                        end do
                    end do
                case (6)
                    !
                    !==> case ityp=6 ,  v odd , w even
                    !
                    !     case m=0
                    !
                    do k=1, nt
                        do i=1, imid
                            do np1=2, ndo2, 2
                                cr(1, np1, k) = cr(1, np1, k)-vb(i, np1)*we(i, 1, k)
                            end do
                        end do
                    end do

                    do k=1, nt
                        do i=1, imm1
                            do np1=3, ndo1, 2
                                br(1, np1, k) = br(1, np1, k)+vb(i, np1)*vo(i, 1, k)
                            end do
                        end do
                    end do
                    !
                    !==> case m = 1 through nlat-1
                    !
                    if (mmax < 2) return

                    do mp1=2, mmax
                        m = mp1-1
                        mb = m*nlat-(m*(m+1))/2
                        mp2 = mp1+1

                        if (mp1 <= ndo1) then

                            do k=1, nt
                                do i=1, imm1
                                    do np1=mp1, ndo1, 2
                                        br(mp1, np1, k) = br(mp1, np1, k)+vb(i, np1+mb)*vo(i, 2*mp1-2, k) &
                                            +wb(i, np1+mb)*we(i, 2*mp1-1, k)
                                        bi(mp1, np1, k) = bi(mp1, np1, k)+vb(i, np1+mb)*vo(i, 2*mp1-1, k) &
                                            -wb(i, np1+mb)*we(i, 2*mp1-2, k)
                                    end do
                                end do
                            end do


                            if (mlat /= 0) then
                                do k=1, nt
                                    do np1=mp1, ndo1, 2
                                        br(mp1, np1, k) = br(mp1, np1, k)+wb(imid, np1+mb)*we(imid, 2*mp1-1, k)
                                        bi(mp1, np1, k) = bi(mp1, np1, k)-wb(imid, np1+mb)*we(imid, 2*mp1-2, k)
                                    end do
                                end do
                            end if

                        end if

                        if (mp2 > ndo2) exit

                        do k=1, nt
                            do i=1, imm1
                                do np1=mp2, ndo2, 2
                                    cr(mp1, np1, k) = cr(mp1, np1, k)-vb(i, np1+mb)*we(i, 2*mp1-2, k) &
                                        +wb(i, np1+mb)*vo(i, 2*mp1-1, k)
                                    ci(mp1, np1, k) = ci(mp1, np1, k)-vb(i, np1+mb)*we(i, 2*mp1-1, k) &
                                        -wb(i, np1+mb)*vo(i, 2*mp1-2, k)
                                end do
                            end do
                        end do

                        if (mlat == 0) exit

                        do k=1, nt
                            do np1=mp2, ndo2, 2
                                cr(mp1, np1, k) = cr(mp1, np1, k)-vb(imid, np1+mb)*we(imid, 2*mp1-2, k)
                                ci(mp1, np1, k) = ci(mp1, np1, k)-vb(imid, np1+mb)*we(imid, 2*mp1-1, k)
                            end do
                        end do
                    end do
                case (7)
                    !
                    !==> case ityp=7   v odd, w even, and cr and ci equal zero
                    !
                    !    case m=0
                    !
                    do k=1, nt
                        do i=1, imm1
                            do np1=3, ndo1, 2
                                br(1, np1, k) = br(1, np1, k)+vb(i, np1)*vo(i, 1, k)
                            end do
                        end do
                    end do
                    !
                    !==> case m = 1 through nlat-1
                    !
                    if (mmax < 2) return

                    do mp1=2, mmax
                        m = mp1-1
                        mb = m*nlat-(m*(m+1))/2
                        mp2 = mp1+1

                        if (mp1 > ndo1) exit

                        do k=1, nt
                            do i=1, imm1
                                do np1=mp1, ndo1, 2
                                    br(mp1, np1, k) = br(mp1, np1, k)+vb(i, np1+mb)*vo(i, 2*mp1-2, k) &
                                        +wb(i, np1+mb)*we(i, 2*mp1-1, k)
                                    bi(mp1, np1, k) = bi(mp1, np1, k)+vb(i, np1+mb)*vo(i, 2*mp1-1, k) &
                                        -wb(i, np1+mb)*we(i, 2*mp1-2, k)
                                end do
                            end do
                        end do


                        if (mlat == 0) exit

                        do  k=1, nt
                            do np1=mp1, ndo1, 2
                                br(mp1, np1, k) = br(mp1, np1, k)+wb(imid, np1+mb)*we(imid, 2*mp1-1, k)
                                bi(mp1, np1, k) = bi(mp1, np1, k)-wb(imid, np1+mb)*we(imid, 2*mp1-2, k)
                            end do
                        end do
                    end do
                case (8)
                    !
                    !==> case ityp=8   v odd, w even, and both br and bi equal zero
                    !
                    !     case m=0
                    !
                    do k=1, nt
                        do i=1, imid
                            do np1=2, ndo2, 2
                                cr(1, np1, k) = cr(1, np1, k)-vb(i, np1)*we(i, 1, k)
                            end do
                        end do
                    end do
                    !
                    !==> case m = 1 through nlat-1
                    !
                    if (mmax < 2) return

                    do mp1=2, mmax
                        m = mp1-1
                        mb = m*nlat-(m*(m+1))/2
                        mp2 = mp1+1

                        if (mp2 > ndo2) exit

                        do k=1, nt
                            do i=1, imm1
                                do np1=mp2, ndo2, 2
                                    cr(mp1, np1, k) = cr(mp1, np1, k)-vb(i, np1+mb)*we(i, 2*mp1-2, k) &
                                        +wb(i, np1+mb)*vo(i, 2*mp1-1, k)
                                    ci(mp1, np1, k) = ci(mp1, np1, k)-vb(i, np1+mb)*we(i, 2*mp1-1, k) &
                                        -wb(i, np1+mb)*vo(i, 2*mp1-2, k)
                                end do
                            end do
                        end do

                        if (mlat == 0) exit

                        do  k=1, nt
                            do  np1=mp2, ndo2, 2
                                cr(mp1, np1, k) = cr(mp1, np1, k)-vb(imid, np1+mb)*we(imid, 2*mp1-2, k)
                                ci(mp1, np1, k) = ci(mp1, np1, k)-vb(imid, np1+mb)*we(imid, 2*mp1-1, k)
                            end do
                        end do
                    end do
            end select

        end subroutine vhags1

    end subroutine vhags

    subroutine vhagsi(nlat, nlon, wvhags, lvhags, dwork, ldwork, ierror)

        integer (ip) :: ierror
        integer (ip) :: imid
        integer (ip) :: iw1
        integer (ip) :: iw2
        integer (ip) :: iw3
        integer (ip) :: iw4
        integer (ip) :: jw1
        integer (ip) :: jw2
        integer (ip) :: jw3
        integer (ip) :: ldwork
        integer (ip) :: lmn
        integer (ip) :: lvhags
        integer (ip) :: nlat
        integer (ip) :: nlon
        real (wp) :: wvhags
        dimension wvhags(lvhags)
        real dwork(ldwork)

        type (HFFTpack) :: hfft

        ierror = 1
        if (nlat < 3) return
        ierror = 2
        if (nlon < 1) return
        ierror = 3
        imid = (nlat+1)/2
        lmn = (nlat*(nlat+1))/2
        if (lvhags < 2*(imid*lmn)+nlon+15) return
        ierror = 4
        !     if (ldwork.lt.nlat*(3*nlat+9)+2) return
        if (ldwork<(nlat*(3*nlat+9)+2)/2) return
        ierror = 0
        jw1 = 1
        jw2 = jw1+imid*lmn
        jw3 = jw2+imid*lmn
        iw1 = 1
        iw2 = iw1+nlat
        iw3 = iw2+nlat
        iw4 = iw3+3*imid*nlat
        call vhgai1(nlat, imid, wvhags(jw1), wvhags(jw2), &
            dwork(iw1), dwork(iw2), dwork(iw3), dwork(iw4))

        call hfft%initialize(nlon, wvhags(jw3))

    contains

        subroutine vhgai1(nlat, imid, vb, wb, dthet, dwts, dpbar, work)
            integer (ip) :: i
            integer (ip) :: id
            integer (ip) :: ierror
            integer (ip) :: imid
            integer (ip) :: ix
            integer (ip) :: iy
            integer (ip) :: lwk
            integer (ip) :: m
            integer (ip) :: mn
            integer (ip) :: n
            integer (ip) :: nlat
            integer (ip) :: nm
            integer (ip) :: np
            integer (ip) :: nz
            real (wp) :: vb
            real (wp) :: wb
            dimension vb(imid, *), wb(imid, *)
            real (wp), parameter :: PI = acos(-1.0_wp)
            real (wp) :: abel, bbel, cbel, dcf
            real (wp) :: dpbar(imid, nlat, 3), dthet(*), dwts(*), work(*)
            real (wp) :: dummy_variable

            type (SpherepackAux) :: sphere_aux

            !     lwk = 4*nlat*(nlat+2)
            lwk = nlat*(nlat+2)

            call gaqd(nlat, dthet, dwts, dummy_variable, lwk, ierror)
            !
            !     compute associated legendre functions
            !
            !     compute m=n=0 legendre polynomials for all theta(i)
            !
            dpbar(1:imid, 1, 1) = cos(PI/4)
            vb(1:imid, 1) = 0.0
            wb(1:imid, 1) = 0.0
            !
            !==> main loop for remaining vb, and wb
            !
            do n=1, nlat-1
                nm = mod(n-2, 3)+1
                nz = mod(n-1, 3)+1
                np = mod(n, 3)+1
                !
                !==> compute dpbar for m=0
                !
                call sphere_aux%dnlfk(0, n, work)
                mn = indx(0, n, nlat)
                do i=1, imid
                    call sphere_aux%dnlft(0, n, dthet(i), work, dpbar(i, 1, np))
                end do
                !
                !==> compute dpbar for m=1
                !
                call sphere_aux%dnlfk(1, n, work)
                mn = indx(1, n, nlat)
                do i=1, imid
                    call sphere_aux%dnlft(1, n, dthet(i), work, dpbar(i, 2, np))
                !      pbar(i, mn) = dpbar(i, 2, np)
                end do

                !
                !==> compute and store dpbar for m=2, n
                !
                if (n >= 2) then
                    do m=2, n
                        abel = sqrt(real((2*n+1)*(m+n-2)*(m+n-3))/ &
                            real((2*n-3)*(m+n-1)*(m+n)))
                        bbel = sqrt(real((2*n+1)*(n-m-1)*(n-m))/ &
                            real((2*n-3)*(m+n-1)*(m+n)))
                        cbel = sqrt(real((n-m+1)*(n-m+2))/ &
                            real((m+n-1)*(m+n)))
                        id = indx(m, n, nlat)

                        if (m >= n-1) then
                            dpbar(1:imid, m+1, np) = &
                                abel*dpbar(1:imid, m-1, nm)-cbel*dpbar(1:imid, m-1, np)
                        else
                            dpbar(1:imid, m+1, np) = &
                                abel*dpbar(1:imid, m-1, nm)+bbel*dpbar(1:imid, m+1, nm) &
                                -cbel*dpbar(1:imid, m-1, np)
                        end if
                    end do
                end if
                !
                !     compute the derivative of the functions
                !
                ix = indx(0, n, nlat)
                iy = indx(n, n, nlat)
                vb(1:imid, ix) = -dpbar(1:imid, 2, np)*dwts(1:imid)
                vb(1:imid, iy) = dpbar(1:imid, n, np)/sqrt(real(2*(n+1)))*dwts(1:imid)

                if (n==1) then
                    !
                    !==> compute the vector harmonic w(theta) = m*pbar/cos(theta)
                    !
                    !     set wb=0 for m=0
                    !
                    ix = indx(0, n, nlat)
                    wb(1:imid, ix) = 0.0
                else
                    dcf = sqrt(real(4*n*(n+1)))
                    do m=1, n-1
                        ix = indx(m, n, nlat)
                        abel = sqrt(real((n+m)*(n-m+1)))/dcf
                        bbel = sqrt(real((n-m)*(n+m+1)))/dcf
                        vb(1:imid, ix) = &
                            (abel*dpbar(1:imid, m, np)-bbel*dpbar(1:imid, m+2, np))&
                            * dwts(1:imid)
                    end do
                end if
                !
                !==> compute wb for m=1, n
                !
                dcf = sqrt(real(n+n+1)/real(4*n*(n+1)*(n+n-1)))
                do m=1, n
                    ix = indx(m, n, nlat)
                    abel = dcf*sqrt(real((n+m)*(n+m-1)))
                    bbel = dcf*sqrt(real((n-m)*(n-m-1)))
                    if (m >= n-1) then
                        wb(1:imid, ix) = abel * dpbar(1:imid, m, nz) * dwts(1:imid)
                    else
                        wb(1:imid, ix) = &
                            (abel*dpbar(1:imid, m, nz) + bbel*dpbar(1:imid, m+2, nz))&
                            * dwts(1:imid)
                    end if
                end do
            end do

        end subroutine vhgai1



        pure function indx(m, n, nlat) result (return_value)
            !----------------------------------------------------------------------
            ! Dictionary: calling arguments
            !----------------------------------------------------------------------
            integer, intent (in) :: m
            integer, intent (in) :: n
            integer, intent (in) :: nlat
            integer              :: return_value
            !----------------------------------------------------------------------

            return_value = m*nlat-(m*(m+1))/2+n+1

        end function indx

    end subroutine vhagsi

end module module_vhags
