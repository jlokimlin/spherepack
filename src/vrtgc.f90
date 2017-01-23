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
! ... file vrtgc.f
!
!     this file includes documentation and code for
!     subroutine divec          i
!
! ... files which must be loaded with vrtgc.f
!
!     type_SpherepackAux.f, type_RealPeriodicTransform.f, vhagc.f, shsgc.f, compute_gaussian_latitudes_and_weights.f
!
!     subroutine vrtgc(nlat, nlon, isym, nt, vort, ivrt, jvrt, cr, ci, mdc, ndc, 
!    +                 wshsgc, lshsgc, work, lwork, ierror)
!
!     given the vector spherical harmonic coefficients cr and ci, precomputed
!     by subroutine vhagc for a vector field (v, w), subroutine vrtgc
!     computes the vorticity of the vector field in the scalar array
!     vort.  vort(i, j) is the vorticity at the gaussian colatitude
!     theta(i) (see nlat as input parameter) and longitude
!     lambda(j) = (j-1)*2*pi/nlon on the sphere.  i.e., 
!
!            vort(i, j) =  [-dv/dlambda + d(sint*w)/dtheta]/sint
!
!     where sint = sin(theta(i)).  w is the east longitudinal and v
!     is the colatitudinal component of the vector field from which
!     cr, ci were precomputed.  required associated legendre polynomials
!     are recomputed rather than stored as they are in subroutine vrtgs.
!
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
!            than 3. the axisymmetric case corresponds to nlon=1.
!            the efficiency of the computation is improved when nlon
!            is a product of small prime numbers.
!
!
!     isym   a parameter which determines whether the vorticity is
!            computed on the full or half sphere as follows:
!
!      = 0
!            the symmetries/antsymmetries described in isym=1, 2 below
!            do not exist in (v, w) about the equator.  in this case the
!            vorticity is neither symmetric nor antisymmetric about
!            the equator.  the vorticity is computed on the entire
!            sphere.  i.e., in the array vort(i, j) for i=1, ..., nlat and
!            j=1, ..., nlon.
!
!      = 1
!            w is antisymmetric and v is symmetric about the equator.
!            in this case the vorticity is symmetyric about the
!            equator and is computed for the northern hemisphere
!            only.  i.e., if nlat is odd the vorticity is computed
!            in the array vort(i, j) for i=1, ..., (nlat+1)/2 and for
!            j=1, ..., nlon.  if nlat is even the vorticity is computed
!            in the array vort(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
!
!      = 2
!            w is symmetric and v is antisymmetric about the equator
!            in this case the vorticity is antisymmetric about the
!            equator and is computed for the northern hemisphere
!            only.  i.e., if nlat is odd the vorticity is computed
!            in the array vort(i, j) for i=1, ..., (nlat+1)/2 and for
!            j=1, ..., nlon.  if nlat is even the vorticity is computed
!            in the array vort(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
!
!
!      nt    nt is the number of scalar and vector fields.  some
!            computational efficiency is obtained for multiple fields.
!            in the program that calls vrtgc, the arrays cr, ci, and vort
!            can be three dimensional corresponding to an indexed multiple
!            vector field.  in this case multiple scalar synthesis will
!            be performed to compute the vorticity for each field.  the
!            third index is the synthesis index which assumes the values
!            k=1, ..., nt.  for a single synthesis set nt = 1.  the
!            description of the remaining parameters is simplified by
!            assuming that nt=1 or that all the arrays are two dimensional.
!
!     ivrt   the first dimension of the array vort as it appears in
!            the program that calls vrtgc. if isym = 0 then ivrt
!            must be at least nlat.  if isym = 1 or 2 and nlat is
!            even then ivrt must be at least nlat/2. if isym = 1 or 2
!            and nlat is odd then ivrt must be at least (nlat+1)/2.
!
!     jvrt   the second dimension of the array vort as it appears in
!            the program that calls vrtgc. jvrt must be at least nlon.
!
!    cr, ci   two or three dimensional arrays (see input parameter nt)
!            that contain vector spherical harmonic coefficients
!            of the vector field (v, w) as computed by subroutine vhagc.
!     ***    cr and ci must be computed by vhagc prior to calling
!            vrtgc.
!
!      mdc   the first dimension of the arrays cr and ci as it
!            appears in the program that calls vrtgc. mdc must be at
!            least min(nlat, nlon/2) if nlon is even or at least
!            min(nlat, (nlon+1)/2) if nlon is odd.
!
!      ndc   the second dimension of the arrays cr and ci as it
!            appears in the program that calls vrtgc. ndc must be at
!            least nlat.
!
!   wshsgc   an array which must be initialized by subroutine shsgci.
!            once initialized, 
!            wshsgc can be used repeatedly by vrtgc as long as nlon
!            and nlat remain unchanged.  wshsgc must not be altered
!            between calls of vrtgc
!
!   lshsgc   the dimension of the array wshsgc   as it appears in the
!            program that calls vrtgc. define
!
!               l1 = min(nlat, (nlon+2)/2) if nlon is even or
!               l1 = min(nlat, (nlon+1)/2) if nlon is odd
!
!            and
!
!               l2 = nlat/2        if nlat is even or
!               l2 = (nlat+1)/2    if nlat is odd
!
!            then lshsgc must be at least
!
!               nlat*(2*l2+3*l1-2)+3*l1*(1-l1)/2+nlon+15
!
!
!     work   a work array that does not have to be saved.
!
!    lwork   the dimension of the array work as it appears in the
!            program that calls vrtgc. define
!
!               l1 = min(nlat, nlon/2) if nlon is even or
!               l1 = min(nlat, (nlon+1)/2) if nlon is odd
!
!            and
!
!               l2 = nlat/2        if nlat is even or
!               l2 = (nlat+1)/2    if nlat is odd.
!
!            if isym is zero then lwork must be at least
!
!               nlat*(nlon*nt+max(3*l2, nlon) + 2*nt*l1+1)
!
!            if isym is not zero then lwork must be at least
!
!               l2*(nlon*nt+max(3*nlat, nlon)) + nlat*(2*nt*l1+1)
!
!
!     **************************************************************
!
!     output parameters
!
!
!     vort   a two or three dimensional array (see input parameter nt)
!            that contains the vorticity of the vector field (v, w)
!            whose coefficients cr, ci where computed by subroutine vhagc.
!            vort(i, j) is the vorticity at the gaussian colatitude point
!            theta(i) and longitude point lambda(j) = (j-1)*2*pi/nlon.
!            the index ranges are defined above at the input parameter
!            isym.
!
!
!   ierror   an error parameter which indicates fatal errors with input
!            parameters when returned positive.
!          = 0  no errors
!          = 1  error in the specification of nlat
!          = 2  error in the specification of nlon
!          = 3  error in the specification of isym
!          = 4  error in the specification of nt
!          = 5  error in the specification of ivrt
!          = 6  error in the specification of jvrt
!          = 7  error in the specification of mdc
!          = 8  error in the specification of ndc
!          = 9  error in the specification of lshsgc
!          = 10 error in the specification of lwork
! **********************************************************************
!                                                                              
!
module module_vrtgc

    use spherepack_precision, only: &
        wp, & ! working precision
        ip ! integer precision

    use scalar_synthesis_routines, only: &
        shsgc

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: vrtgc

contains

    subroutine vrtgc(nlat, nlon, isym, nt, vort, ivrt, jvrt, cr, ci, mdc, ndc, &
        wshsgc, lshsgc, work, lwork, ierror)

        real(wp) :: ci
        real(wp) :: cr
        integer(ip) :: ia
        integer(ip) :: ib
        integer(ip) :: ierror
        integer(ip) :: imid
        integer(ip) :: is
        integer(ip) :: isym
        integer(ip) :: ivrt
        integer(ip) :: iwk
        integer(ip) :: jvrt
        integer(ip) :: l1
        integer(ip) :: l2
        integer(ip) :: lpimn
        integer(ip) :: ls
        integer(ip) :: lshsgc
        integer(ip) :: lwk
        integer(ip) :: lwkmin
        integer(ip) :: lwmin
        integer(ip) :: lwork
        integer(ip) :: mab
        integer(ip) :: mdc
        integer(ip) :: mmax
        integer(ip) :: mn
        integer(ip) :: ndc
        integer(ip) :: nlat
        integer(ip) :: nln
        integer(ip) :: nlon
        integer(ip) :: nt
        real(wp) :: vort
        real(wp) :: work
        real(wp) :: wshsgc

        dimension vort(ivrt, jvrt, nt), cr(mdc, ndc, nt), ci(mdc, ndc, nt)
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
        if ((isym == 0 .and. ivrt<nlat) .or. &
            (isym>0 .and. ivrt<imid)) return
        ierror = 6
        if (jvrt < nlon) return
        ierror = 7
        if (mdc < min(nlat, (nlon+1)/2)) return
        mmax = min(nlat, (nlon+2)/2)
        ierror = 8
        if (ndc < nlat) return
        ierror = 9
        imid = (nlat+1)/2
        lpimn = (imid*mmax*(nlat+nlat-mmax+1))/2
        l1 = min(nlat, (nlon+2)/2)
        l2 = (nlat+1)/2
        lwmin = nlat*(2*l2+3*l1-2)+3*l1*(1-l1)/2+nlon+15
        if (lshsgc < lwmin) return
        ierror = 10
        !
        !     verify unsaved work space (add to what shses requires, file f3)
        !
        !
        !     set first dimension for a, b (as requried by shses)
        !
        mab = min(nlat, nlon/2+1)
        mn = mab*nlat*nt
        ls = nlat
        if (isym > 0) ls = imid
        nln = nt*ls*nlon
        !     if (lwork.lt. nln+ls*nlon+2*mn+nlat) return
        l1 = min(nlat, (nlon+2)/2)
        l2 = (nlat+1)/2
        if (isym == 0) then
            lwkmin =  nlat*(nt*nlon+max(3*l2, nlon)+2*nt*l1+1)
        else
            lwkmin = l2*(nt*nlon+max(3*nlat, nlon)) + nlat*(2*nt*l1+1)
        end if
        if (lwork < lwkmin) return
        ierror = 0
        !
        !     set work space pointers
        !
        ia = 1
        ib = ia+mn
        is = ib+mn
        iwk = is+nlat
        lwk = lwork-2*mn-nlat
        call vrtgc1(nlat, nlon, isym, nt, vort, ivrt, jvrt, cr, ci, mdc, ndc, &
            work(ia), work(ib), mab, work(is), wshsgc, lshsgc, work(iwk), lwk, &
            ierror)

    contains

        subroutine vrtgc1(nlat, nlon, isym, nt, vort, ivrt, jvrt, cr, ci, mdc, ndc, &
            a, b, mab, sqnn, wsav, lwsav, wk, lwk, ierror)

            real(wp) :: a
            real(wp) :: b
            real(wp) :: ci
            real(wp) :: cr
            real(wp) :: fn
            integer(ip) :: ierror
            integer(ip) :: isym
            integer(ip) :: ivrt
            integer(ip) :: jvrt
            integer(ip) :: k
            integer(ip) :: lwk
            integer(ip) :: lwsav
            integer(ip) :: m
            integer(ip) :: mab
            integer(ip) :: mdc
            integer(ip) :: mmax
            integer(ip) :: n
            integer(ip) :: ndc
            integer(ip) :: nlat
            integer(ip) :: nlon
            integer(ip) :: nt
            real(wp) :: sqnn
            real(wp) :: vort
            real(wp) :: wk
            real(wp) :: wsav
            dimension vort(ivrt, jvrt, nt), cr(mdc, ndc, nt), ci(mdc, ndc, nt)
            dimension a(mab, nlat, nt), b(mab, nlat, nt), sqnn(nlat)
            dimension wsav(lwsav), wk(lwk)
            !
            !     set coefficient multiplyers
            !
            do n=2, nlat
                fn = real(n - 1)
                sqnn(n) = sqrt(fn * (fn + 1.0))
            end do
            !
            !     compute divergence scalar coefficients for each vector field
            !
            do k=1, nt
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
                    a(1, n, k) = sqnn(n)*cr(1, n, k)
                    b(1, n, k) = sqnn(n)*ci(1, n, k)
                end do
                !
                !     compute m>0 coefficients
                !
                mmax = min(nlat, (nlon+1)/2)
                do m=2, mmax
                    do n=m, nlat
                        a(m, n, k) = sqnn(n)*cr(m, n, k)
                        b(m, n, k) = sqnn(n)*ci(m, n, k)
                    end do
                end do
            end do
            !
            !     synthesize a, b into vort
            !
            call shsgc(nlat, nlon, isym, nt, vort, ivrt, jvrt, a, b, &
                mab, nlat, wsav, lwsav, wk, lwk, ierror)

        end subroutine vrtgc1

    end subroutine vrtgc


end module module_vrtgc
