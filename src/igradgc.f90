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
! ... file igradgc.f
!
!     this file includes documentation and code for
!     subroutine igradgc         i
!
! ... files which must be loaded with igradgc.f
!
!     type_SpherepackAux.f, type_RealPeriodicTransform.f, shsgc.f, vhagc.f
!
!     subroutine igradgc(nlat, nlon, isym, nt, sf, isf, jsf, br, bi, mdb, ndb, 
!    +                   wshsgc, lshsgc, work, lwork, ierror)
!
!     let br, bi, cr, ci be the vector spherical harmonic coefficients
!     precomputed by vhagc for a vector field (v, w).  let (v', w') be
!     the irrotational component of (v, w) (i.e., (v', w') is generated
!     by assuming cr, ci are zero and synthesizing br, bi with vhsgs).
!     then subroutine igradgc computes a scalar field sf such that
!
!            gradient(sf) = (v', w').
!
!     i.e., 
!
!            v'(i, j) = d(sf(i, j))/dtheta          (colatitudinal component of
!                                                 the gradient)
!     and
!
!            w'(i, j) = 1/sint*d(sf(i, j))/dlambda  (east longitudinal component
!                                                 of the gradient)
!
!     at the gaussian colatitude theta(i) (see nlat as input parameter)
!     and longitude lambda(j) = (j-1)*2*pi/nlon where sint = sin(theta(i)).
!
!     note:  for an irrotational vector field (v, w), subroutine igradgc
!     computes a scalar field whose gradient is (v, w).  in ay case, 
!     subroutine igradgc "inverts" the gradient subroutine gradgc.
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
!            nlon = 72 for a five degree grid. nlon must be greater than
!            3.  the efficiency of the computation is improved when nlon
!            is a product of small prime numbers.
!
!
!     isym   a parameter which determines whether the scalar field sf is
!            computed on the full or half sphere as follows:
!
!      = 0
!
!            the symmetries/antsymmetries described in isym=1, 2 below
!            do not exist in (v, w) about the equator.  in this case sf
!            is neither symmetric nor antisymmetric about the equator.
!            sf is computed on the entire sphere.  i.e., in the array
!            sf(i, j) for i=1, ..., nlat and  j=1, ..., nlon
!
!      = 1
!
!            w is antisymmetric and v is symmetric about the equator.
!            in this case sf is antisymmetyric about the equator and
!            is computed for the northern hemisphere only.  i.e., 
!            if nlat is odd sf is computed in the array sf(i, j) for
!            i=1, ..., (nlat+1)/2 and for j=1, ..., nlon.  if nlat is even
!            sf is computed in the array sf(i, j) for i=1, ..., nlat/2
!            and j=1, ..., nlon.
!
!      = 2
!
!            w is symmetric and v is antisymmetric about the equator.
!            in this case sf is symmetyric about the equator and
!            is computed for the northern hemisphere only.  i.e., 
!            if nlat is odd sf is computed in the array sf(i, j) for
!            i=1, ..., (nlat+1)/2 and for j=1, ..., nlon.  if nlat is even
!            sf is computed in the array sf(i, j) for i=1, ..., nlat/2
!            and j=1, ..., nlon.
!
!
!     nt     nt is the number of scalar and vector fields.  some
!            computational efficiency is obtained for multiple fields.
!            the arrays br, bi, and sf can be three dimensional corresponding
!            to an indexed multiple vector field (v, w).  in this case, 
!            multiple scalar synthesis will be performed to compute each
!            scalar field.  the third index for br, bi, and sf is the synthesis
!            index which assumes the values k = 1, ..., nt.  for a single
!            synthesis set nt = 1.  the description of the remaining
!            parameters is simplified by assuming that nt=1 or that br, bi, 
!            and sf are two dimensional arrays.
!
!     isf    the first dimension of the array sf as it appears in
!            the program that calls igradgc. if isym = 0 then isf
!            must be at least nlat.  if isym = 1 or 2 and nlat is
!            even then isf must be at least nlat/2. if isym = 1 or 2
!            and nlat is odd then isf must be at least (nlat+1)/2.
!
!     jsf    the second dimension of the array sf as it appears in
!            the program that calls igradgc. jsf must be at least nlon.
!
!     br, bi  two or three dimensional arrays (see input parameter nt)
!            that contain vector spherical harmonic coefficients
!            of the vector field (v, w) as computed by subroutine vhagc.
!     ***    br, bi must be computed by vhagc prior to calling igradgc.
!
!     mdb    the first dimension of the arrays br and bi as it appears in
!            the program that calls igradgc (and vhagc). mdb must be at
!            least min(nlat, nlon/2) if nlon is even or at least
!            min(nlat, (nlon+1)/2) if nlon is odd.
!
!     ndb    the second dimension of the arrays br and bi as it appears in
!            the program that calls igradgc (and vhagc). ndb must be at
!            least nlat.
!
!
!  wshsgc    an array which must be initialized by subroutine shsgci.
!            once initialized, 
!            wshsgc can be used repeatedly by igradgc as long as nlon
!            and nlat remain unchanged.  wshsgc must not be altered
!            between calls of igradgc.
!
!
!  lshsgc    the dimension of the array wshsgc as it appears in the
!            program that calls igradgc. define
!
!               l1 = min(nlat, (nlon+2)/2) if nlon is even or
!               l1 = min(nlat, (nlon+1)/2) if nlon is odd
!
!            and
!
!               l2 = nlat/2        if nlat is even or
!               l2 = (nlat+1)/2    if nlat is odd.
!
!
!            then lshsgc must be at least
!
!               nlat*(2*l2+3*l1-2)+3*l1*(1-l1)/2+nlon+15
!
!
!     work   a work array that does not have to be saved.
!
!     lwork  the dimension of the array work as it appears in the
!            program that calls igradgc  define
!
!               l2 = nlat/2                    if nlat is even or
!               l2 = (nlat+1)/2                if nlat is odd
!               l1 = min(nlat, (nlon+2)/2) if nlon is even or
!               l1 = min(nlat, (nlon+1)/2) if nlon is odd

!            if isym is zero then lwork must be at least
!
!               nlat*(nlon*nt+max(3*l2, nlon)+2*nt*l1+1)
!
!            if isym is not zero then lwork must be at least
!
!               l2*(nlon*nt+max(3*nlat, nlon)) + nlat*(2*nt*l1+1)
!
!
!
!     **************************************************************
!
!     output parameters
!
!
!     sf    a two or three dimensional array (see input parameter nt) that
!           contain a scalar field whose gradient is the irrotational
!           component of the vector field (v, w).  the vector spherical
!           harmonic coefficients br, bi were precomputed by subroutine
!           vhagc.  sf(i, j) is given at the gaussian colatitude theta(i)
!           and longitude lambda(j) = (j-1)*2*pi/nlon.  the index ranges
!           are defined at input parameter isym.
!
!
!  ierror   = 0  no errors
!           = 1  error in the specification of nlat
!           = 2  error in the specification of nlon
!           = 3  error in the specification of isym
!           = 4  error in the specification of nt
!           = 5  error in the specification of isf
!           = 6  error in the specification of jsf
!           = 7  error in the specification of mdb
!           = 8  error in the specification of ndb
!           = 9  error in the specification of lshsgc
!           = 10 error in the specification of lwork
!
! **********************************************************************
!
module module_igradgc

    use spherepack_precision, only: &
        wp, & ! working precision
        ip ! integer precision

    use scalar_synthesis_routines, only: &
        shsgc

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    public :: igradgc

contains

    subroutine igradgc(nlat, nlon, isym, nt, sf, isf, jsf, br, bi, mdb, ndb, &
        wshsgc, lshsgc, work, lwork, ierror)

        real(wp) :: bi
        real(wp) :: br
        integer(ip) :: ia
        integer(ip) :: ib
        integer(ip) :: ierror
        integer(ip) :: imid
        integer(ip) :: is
        integer(ip) :: isf
        integer(ip) :: isym
        integer(ip) :: iwk
        integer(ip) :: jsf
        integer(ip) :: l1
        integer(ip) :: l2
        integer(ip) :: liwk
        integer(ip) :: ls
        integer(ip) :: lshsgc
        integer(ip) :: lwkmin
        integer(ip) :: lwork
        integer(ip) :: mab
        integer(ip) :: mdb
        integer(ip) :: mmax
        integer(ip) :: mn
        integer(ip) :: ndb
        integer(ip) :: nlat
        integer(ip) :: nln
        integer(ip) :: nlon
        integer(ip) :: nt
        real(wp) :: sf
        real(wp) :: work
        real(wp) :: wshsgc
        dimension sf(isf, jsf, nt)
        dimension br(mdb, ndb, nt), bi(mdb, ndb, nt)
        dimension wshsgc(lshsgc), work(lwork)
        !
        ! Check input arguments
        !
        ierror = 1
        if (nlat < 3) return
        ierror = 2
        if (nlon < 4) return
        ierror = 3
        if (isym < 0 .or. isym > 2) return
        ierror = 4
        if (nt < 0) return
        ierror = 5
        imid = (nlat+1)/2
        if ((isym == 0 .and. isf<nlat) .or. &
            (isym /= 0 .and. isf<imid)) return
        ierror = 6
        if (jsf < nlon) return
        ierror = 7
        mmax = min(nlat, (nlon+2)/2)
        if (mdb < min(nlat, (nlon+1)/2)) return
        ierror = 8
        if (ndb < nlat) return
        ierror = 9
        !
        !     verify saved work space length
        !
        l2 = (nlat+mod(nlat, 2))/2
        l1 = min((nlon+2)/2, nlat)
        if (lshsgc < nlat*(2*l2+3*l1-2)+3*l1*(1-l1)/2+nlon+15)return
        ierror = 10
        !
        !     set minimum and verify unsaved work space
        !
        ls = nlat
        if (isym > 0) ls = imid
        nln = nt*ls*nlon
        !
        !     set first dimension for a, b (as requried by shsgc)
        !
        mab = min(nlat, nlon/2+1)
        mn = mab*nlat*nt
        !     lwkmin = nln+ls*nlon+2*mn+nlat
        if (isym == 0) then
            lwkmin = nlat*(nt*nlon+max(3*l2, nlon)+2*nt*l1+1)
        else
            lwkmin = l2*(nt*nlon+max(3*nlat, nlon))+nlat*(2*nt*l1+1)
        end if
        if (lwork < lwkmin) return
        ierror = 0
        !
        !     set work space pointers
        !
        ia = 1
        ib = ia + mn
        is = ib + mn
        iwk = is + nlat
        liwk = lwork-2*mn-nlat

        call igrdgc1(nlat, nlon, isym, nt, sf, isf, jsf, work(ia), work(ib), mab, &
            work(is), mdb, ndb, br, bi, wshsgc, lshsgc, work(iwk), liwk, ierror)


    contains


        subroutine igrdgc1(nlat, nlon, isym, nt, sf, isf, jsf, a, b, mab, &
            sqnn, mdb, ndb, br, bi, wsav, lsav, wk, lwk, ierror)

            real(wp) :: a
            real(wp) :: b
            real(wp) :: bi
            real(wp) :: br
            real(wp) :: fn
            integer(ip) :: ierror
            integer(ip) :: isf
            integer(ip) :: isym
            integer(ip) :: jsf
            integer(ip) :: k
            integer(ip) :: lsav
            integer(ip) :: lwk
            integer(ip) :: m
            integer(ip) :: mab
            integer(ip) :: mdb
            integer(ip) :: mmax
            integer(ip) :: n
            integer(ip) :: ndb
            integer(ip) :: nlat
            integer(ip) :: nlon
            integer(ip) :: nt
            real(wp) :: sf
            real(wp) :: sqnn
            real(wp) :: wk
            real(wp) :: wsav
            dimension sf(isf, jsf, nt)
            dimension br(mdb, ndb, nt), bi(mdb, ndb, nt), sqnn(nlat)
            dimension a(mab, nlat, nt), b(mab, nlat, nt)
            dimension wsav(lsav), wk(lwk)
            !
            ! Preset coefficient multiplyers in vector
            !
            do n=2, nlat
                fn = real(n - 1)
                sqnn(n) = 1.0/sqrt(fn * (fn + 1.0))
            end do
            !
            !     set upper limit for vector m subscript
            !
            mmax = min(nlat, (nlon+1)/2)
            !
            !     compute multiple scalar field coefficients
            !
            do k=1, nt
                !
                !     preset to 0.0
                !
                do n=1, nlat
                    do m=1, mab
                        a(m, n, k) = 0.0
                        b(m, n, k) = 0.0
                    end do
                end do
                !
                !     compute m=0 coefficients
                !
                do n=2, nlat
                    a(1, n, k) = br(1, n, k)*sqnn(n)
                    b(1, n, k)= bi(1, n, k)*sqnn(n)
                end do
                !
                !     compute m>0 coefficients
                !
                do m=2, mmax
                    do n=m, nlat
                        a(m, n, k) = sqnn(n)*br(m, n, k)
                        b(m, n, k) = sqnn(n)*bi(m, n, k)
                    end do
                end do
            end do
            !
            !     scalar sythesize a, b into sf
            !
            call shsgc(nlat, nlon, isym, nt, sf, isf, jsf, a, b, mab, nlat, wsav, &
                lsav, wk, lwk, ierror)

        end subroutine igrdgc1

    end subroutine igradgc

end module module_igradgc
