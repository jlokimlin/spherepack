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
! ... file divgc.f
!
!     this file includes documentation and code for
!     subroutine divgc          i
!
! ... files which must be loaded with divgc.f
!
!     type_SpherepackAux.f, type_HFFTpack.f, vhagc.f, shsgc.f, compute_gaussian_latitudes_and_weights.f
!
!
!     subroutine divgc(nlat, nlon, isym, nt, dv, idv, jdv, br, bi, mdb, ndb, 
!    +                 wshsgc, lshsgc, work, lwork, ierror)
!
!     given the vector spherical harmonic coefficients br and bi, precomputed
!     by subroutine vhagc for a vector field (v, w), subroutine divgc
!     computes the divergence of the vector field in the scalar array dv.
!     dv(i, j) is the divergence at the gaussian colatitude point theta(i)
!     (see nlat as input parameter) and east longitude
!
!            lambda(j) = (j-1)*2*pi/nlon
!
!     on the sphere.  i.e.
!
!            dv(i, j) = 1/sint*[ d(sint*v(i, j))/dtheta + d(w(i, j))/dlambda ]
!
!     where sint = sin(theta(i)).  w is the east longitudinal and v
!     is the colatitudinal component of the vector field from which
!     br, bi were precomputed.  required associated legendre polynomials
!     are recomputed rather than stored as they are in subroutine divgs.
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
!            than zero. the axisymmetric case corresponds to nlon=1.
!            the efficiency of the computation is improved when nlon
!            is a product of small prime numbers.
!
!
!     isym   a parameter which determines whether the divergence is
!            computed on the full or half sphere as follows:
!
!      = 0
!
!            the symmetries/antsymmetries described in isym=1, 2 below
!            do not exist in (v, w) about the equator.  in this case the
!            divergence is neither symmetric nor antisymmetric about
!            the equator.  the divergence is computed on the entire
!            sphere.  i.e., in the array dv(i, j) for i=1, ..., nlat and
!            j=1, ..., nlon.
!
!      = 1
!
!            w is antisymmetric and v is symmetric about the equator.
!            in this case the divergence is antisymmetyric about
!            the equator and is computed for the northern hemisphere
!            only.  i.e., if nlat is odd the divergence is computed
!            in the array dv(i, j) for i=1, ..., (nlat+1)/2 and for
!            j=1, ..., nlon.  if nlat is even the divergence is computed
!            in the array dv(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
!
!      = 2
!            w is symmetric and v is antisymmetric about the equator
!            in this case the divergence is symmetyric about the
!            equator and is computed for the northern hemisphere
!            only.  i.e., if nlat is odd the divergence is computed
!            in the array dv(i, j) for i=1, ..., (nlat+1)/2 and for
!            j=1, ..., nlon.  if nlat is even the divergence is computed
!            in the array dv(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
!
!
!     nt     nt is the number of scalar and vector fields.  some
!            computational efficiency is obtained for multiple fields.
!            in the program that calls divgc, the arrays br, bi, and dv
!            can be three dimensional corresponding to an indexed multiple
!            vector field.  in this case multiple scalar synthesis will
!            be performed to compute the divergence for each field.  the
!            third index is the synthesis index which assumes the values
!            k=1, ..., nt.  for a single synthesis set nt = 1.  the
!            description of the remaining parameters is simplified by
!            assuming that nt=1 or that all the arrays are two dimensional.
!
!     idv    the first dimension of the array dv as it appears in
!            the program that calls divgc. if isym = 0 then idv
!            must be at least nlat.  if isym = 1 or 2 and nlat is
!            even then idv must be at least nlat/2. if isym = 1 or 2
!            and nlat is odd then idv must be at least (nlat+1)/2.
!
!     jdv    the second dimension of the array dv as it appears in
!            the program that calls divgc. jdv must be at least nlon.
!
!     br, bi  two or three dimensional arrays (see input parameter nt)
!            that contain vector spherical harmonic coefficients
!            of the vector field (v, w) as computed by subroutine vhagc.
!     ***    br and bi must be computed by vhagc prior to calling
!            divgc.
!
!     mdb    the first dimension of the arrays br and bi as it
!            appears in the program that calls divgc. mdb must be at
!            least min(nlat, nlon/2) if nlon is even or at least
!            min(nlat, (nlon+1)/2) if nlon is odd.
!
!     ndb    the second dimension of the arrays br and bi as it
!            appears in the program that calls divgc. ndb must be at
!            least nlat.
!
!
!     wshsgc an array which must be initialized by subroutine shsgci
!            once initialized, wshsgc can be used repeatedly by divgc
!            as long as nlon and nlat remain unchanged.  wshsgc must
!            not be altered between calls of divgc.
!
!
!     lshsgc the dimension of the array wshsgc as it appears in the
!            program that calls divgc. define
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
!     work   a work array that does not have to be saved.
!
!     lwork  the dimension of the array work as it appears in the
!            program that calls divgc. define
!
!               l1 = min(nlat, (nlon+2)/2) if nlon is even or
!               l1 = min(nlat, (nlon+1)/2) if nlon is odd
!
!            and
!
!               l2 = nlat/2                    if nlat is even or
!               l2 = (nlat+1)/2                if nlat is odd
!
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
!    dv     a two or three dimensional array (see input parameter nt)
!           that contains the divergence of the vector field (v, w)
!           whose coefficients br, bi where computed by subroutine
!           vhagc.  dv(i, j) is the divergence at the gaussian colatitude
!           point theta(i) and longitude point lambda(j) = (j-1)*2*pi/nlon.
!           the index ranges are defined above at the input parameter
!           isym.
!
!
!    ierror = 0  no errors
!           = 1  error in the specification of nlat
!           = 2  error in the specification of nlon
!           = 3  error in the specification of isym
!           = 4  error in the specification of nt
!           = 5  error in the specification of idv
!           = 6  error in the specification of jdv
!           = 7  error in the specification of mdb
!           = 8  error in the specification of ndb
!           = 9  error in the specification of lshsgc
!           = 10 error in the specification of lwork
!
!
submodule(divergence_routines) divergence_gaussian_grid

contains

    module subroutine divgc(nlat, nlon, isym, nt, dv, idv, jdv, br, bi, mdb, ndb, &
        wshsgc, lshsgc, work, lwork, ierror)
        real(wp) :: dv(idv, jdv, nt), br(mdb, ndb, nt), bi(mdb, ndb, nt)
        integer(ip) :: ia
        integer(ip) :: ib
        integer(ip) :: idv
        integer(ip) :: ierror
        integer(ip) :: imid
        integer(ip) :: is
        integer(ip) :: isym
        integer(ip) :: iwk
        integer(ip) :: jdv
        integer(ip) :: l1
        integer(ip) :: l2
        integer(ip) :: lpimn
        integer(ip) :: ls
        integer(ip) :: lshsgc
        integer(ip) :: lwk
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
        real(wp) :: wshsgc(lshsgc), work(lwork)
        !
        !     check input parameters
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
        if ((isym == 0 .and. idv<nlat) .or. &
            (isym>0 .and. idv<imid)) return
        ierror = 6
        if (jdv < nlon) return
        ierror = 7
        if (mdb < min(nlat, (nlon+1)/2)) return
        mmax = min(nlat, (nlon+2)/2)
        ierror = 8
        if (ndb < nlat) return
        ierror = 9
        imid = (nlat+1)/2
        lpimn = (imid*mmax*(nlat+nlat-mmax+1))/2
        !     check permanent work space length
        l2 = (nlat+1)/2
        l1 = min((nlon+2)/2, nlat)
        if (lshsgc < nlat*(2*l2+3*l1-2)+3*l1*(1-l1)/2+nlon+15)return
        ierror = 10
        !
        !     verify unsaved work space (add to what shsgc requires)
        !
        ls = nlat
        if (isym > 0) ls = imid
        nln = nt*ls*nlon
        !
        !     set first dimension for a, b (as requried by shsgc)
        !
        mab = min(nlat, nlon/2+1)
        mn = mab*nlat*nt
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
        call divgc_lower_routine(nlat, nlon, isym, nt, dv, idv, jdv, br, bi, mdb, ndb, &
            work(ia), work(ib), mab, work(is), wshsgc, lshsgc, work(iwk), lwk, &
            ierror)

    end subroutine divgc

    subroutine divgc_lower_routine(nlat, nlon, isym, nt, dv, idv, jdv, br, bi, mdb, ndb, &
        a, b, mab, sqnn, wshsgc, lshsgc, wk, lwk, ierror)
        real(wp) :: a
        real(wp) :: b
        real(wp) :: bi
        real(wp) :: br
        real(wp) :: dv
        real(wp) :: fn
        integer(ip) :: idv
        integer(ip) :: ierror
        integer(ip) :: isym
        integer(ip) :: jdv
        integer(ip) :: k
        integer(ip) :: lshsgc
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
        real(wp) :: sqnn
        real(wp) :: wk
        real(wp) :: wshsgc
        dimension dv(idv, jdv, nt), br(mdb, ndb, nt), bi(mdb, ndb, nt)
        dimension a(mab, nlat, nt), b(mab, nlat, nt), sqnn(nlat)
        dimension wshsgc(lshsgc), wk(lwk)
        !
        !     set coefficient multiplyers
        !
        do  n=2, nlat
            fn = real(n - 1, kind=wp)
            sqnn(n) = sqrt(fn * (fn + ONE))
        end do
        !
        !     compute divergence scalar coefficients for each vector field
        !
        do  k=1, nt
            a(1: mab, 1: nlat, k) = ZERO
            b(1: mab, 1: nlat, k) = ZERO
            !
            !     compute m=0 coefficients
            !
            do  n=2, nlat
                a(1, n, k) = -sqnn(n)*br(1, n, k)
                b(1, n, k) = -sqnn(n)*bi(1, n, k)
            end do
            !
            !     compute m>0 coefficients using vector spherepack value for mmax
            !
            mmax = min(nlat, (nlon+1)/2)
            do  m=2, mmax
                do  n=m, nlat
                    a(m, n, k) = -sqnn(n)*br(m, n, k)
                    b(m, n, k) = -sqnn(n)*bi(m, n, k)
                end do
            end do
        end do
        !
        !     synthesize a, b into dv
        !
        call shsgc(nlat, nlon, isym, nt, dv, idv, jdv, a, b, &
            mab, nlat, wshsgc, lshsgc, wk, lwk, ierror)

    end subroutine divgc_lower_routine

end submodule divergence_gaussian_grid
