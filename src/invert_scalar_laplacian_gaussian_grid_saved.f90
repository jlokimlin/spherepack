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
! ... file islapgs.f
!
!     this file includes documentation and code for
!     subroutine islapgs         i
!
! ... files which must be loaded with islapec.f
!
!     type_SpherepackAux.f, type_RealPeriodicTransform.f, shags.f, shsgs.f
!
!     subroutine islapgs(nlat, nlon, isym, nt, xlmbda, sf, ids, jds, a, b, 
!    +mdab, ndab, wshsgs, lshsgs, work, lwork, pertrb, ierror)
!
!     islapgs inverts the laplace or helmholz operator on a Gaussian grid.
!     Given the spherical harmonic coefficients a(m, n) and b(m, n) of the
!     right hand side slap(i, j), islapgc computes a solution sf(i, j) to
!     the following helmhotz equation :
!
!           2                2
!     [d(sf(i, j))/dlambda /sint + d(sint*d(sf(i, j))/dtheta)/dtheta]/sint
!
!                   - xlmbda * sf(i, j) = slap(i, j)
!
!      where sf(i, j) is computed at the Gaussian colatitude point theta(i)
!      (see nlat as an input argument) and longitude
!
!                 lambda(j) = (j-1)*2*pi/nlon
!
!            for i=1, ..., nlat and j=1, ..., nlon.
!
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
!     nlon   the number of distinct longitude points.  nlon determines
!            the grid increment in longitude as 2*pi/nlon. for example
!            nlon = 72 for a five degree grid. nlon must be greater
!            than zero. the axisymmetric case corresponds to nlon=1.
!            the efficiency of the computation is improved when nlon
!            is a product of small prime numbers.
!
!     isym   this parameter should have the same value input to subroutine
!            shags to compute the coefficients a and b for the scalar field
!            slap.  isym is set as follows:
!
!            = 0  no symmetries exist in slap about the equator. scalar
!                 synthesis is used to compute sf on the entire sphere.
!                 i.e., in the array sf(i, j) for i=1, ..., nlat and
!                 j=1, ..., nlon.
!
!           = 1  sf and slap are antisymmetric about the equator. the
!                synthesis used to compute sf is performed on the
!                northern hemisphere only.  if nlat is odd, sf(i, j) is
!                computed for i=1, ..., (nlat+1)/2 and j=1, ..., nlon.  if
!                nlat is even, sf(i, j) is computed for i=1, ..., nlat/2
!                and j=1, ..., nlon.
!
!
!           = 2  sf and slap are symmetric about the equator. the
!                synthesis used to compute sf is performed on the
!                northern hemisphere only.  if nlat is odd, sf(i, j) is
!                computed for i=1, ..., (nlat+1)/2 and j=1, ..., nlon.  if
!                nlat is even, sf(i, j) is computed for i=1, ..., nlat/2
!                and j=1, ..., nlon.
!
!
!     nt     the number of analyses.  in the program that calls islapgs
!            the arrays sf, a, and b can be three dimensional in which
!            case multiple synthesis will be performed.  the third index
!            is the synthesis index which assumes the values k=1, ..., nt.
!            k is also the index for the perturbation array pertrb.
!            for a single analysis set nt=1. the description of the
!            remaining parameters is simplified by assuming that nt=1
!            or that sf, a, b are two dimensional and pertrb is a constant.
!
!   xlmbda   a one dimensional array with nt elements. if xlmbda is
!            is identically zero islapgc solves poisson's equation.
!            if xlmbda > 0.0 islapgc solves the helmholtz equation.
!            if xlmbda < 0.0 the nonfatal error flag ierror=-1 is
!            returned. negative xlambda could result in a division
!            by zero.
!
!   ids      the first dimension of the array sf as it appears in the
!            program that calls islapgs.  if isym = 0 then ids must be at
!            least nlat.  if isym > 0 and nlat is even then ids must be
!            at least nlat/2. if isym > 0 and nlat is odd then ids must
!            be at least (nlat+1)/2.
!
!   jds      the second dimension of the array sf as it appears in the
!            program that calls islapgs. jds must be at least nlon.
!
!
!   a, b      two or three dimensional arrays (see input parameter nt)
!            that contain scalar spherical harmonic coefficients
!            of the scalar field slap as computed by subroutine shags.
!     ***    a, b must be computed by shags prior to calling islapgs.
!
!
!    mdab    the first dimension of the arrays a and b as it appears
!            in the program that calls islapgs.  mdab must be at
!            least min(nlat, (nlon+2)/2) if nlon is even or at least
!            min(nlat, (nlon+1)/2) if nlon is odd.
!
!    ndab    the second dimension of the arrays a and b as it appears
!            in the program that calls islapgs. ndbc must be at least
!            least nlat.
!
!            mdab, ndab should have the same values input to shags to
!            compute the coefficients a and b.
!
!
!    wshsgs  an array which must be initialized by subroutine islapgsi
!            (or equivalently by shsesi).  once initialized, wshsgs
!            can be used repeatedly by islapgs as long as nlat and nlon
!            remain unchanged.  wshsgs  must not be altered between calls
!            of islapgs.
!
!    lshsgs  the dimension of the array wshsgs as it appears in the
!            program that calls islapgs.  let
!
!               l1 = min(nlat, (nlon+2)/2) if nlon is even or
!               l1 = min(nlat, (nlon+1)/2) if nlon is odd
!
!            and
!
!               l2 = nlat/2        if nlat is even or
!               l2 = (nlat+1)/2    if nlat is odd
!
!            then lshsgs must be at least
!
!            nlat*(3*(l1+l2)-2)+(l1-1)*(l2*(2*nlat-l1)-3*l1)/2+nlon+15
!
!     work   a work array that does not have to be saved.
!
!     lwork  the dimension of the array work as it appears in the
!            program that calls islapgs. define
!
!               l2 = nlat/2                    if nlat is even or
!               l2 = (nlat+1)/2                if nlat is odd
!               l1 = min(nlat, (nlon+2)/2) if nlon is even or
!               l1 = min(nlat, (nlon+1)/2) if nlon is odd
!
!            if isym is zero then lwork must be at least
!
!               (nt+1)*nlat*nlon + nlat*(2*nt*l1+1)
!
!            if isym is nonzero lwork must be at least
!
!               (nt+1)*l2*nlon + nlat*(2*nt*l1+1)
!
!
!     **************************************************************
!
!     output parameters
!
!
!    sf      a two or three dimensional arrays (see input parameter nt) that
!            inverts the scalar laplacian in slap.  sf(i, j) is given at
!            the colatitude
!
!                 theta(i) = (i-1)*pi/(nlat-1)
!
!            and longitude
!
!                 lambda(j) = (j-1)*2*pi/nlon
!
!            for i=1, ..., nlat and j=1, ..., nlon.
!
!   pertrb  a one dimensional array with nt elements (see input 
!           parameter nt). in the discription that follows we assume 
!           that nt=1. if xlmbda > 0.0 then pertrb=0.0 is always 
!           returned because the helmholtz operator is invertible.
!           if xlmbda = 0.0 then a solution exists only if a(1, 1)
!           is zero. islapec sets a(1, 1) to zero. the resulting
!           solution sf(i, j) solves poisson's equation with
!           pertrb = a(1, 1)/(2.*sqrt(2.)) subtracted from the
!           right side slap(i, j).
!
!  ierror    a parameter which flags errors in input parameters as follows:
!
!            = 0  no errors detected
!
!            = 1  error in the specification of nlat
!
!            = 2  error in the specification of nlon
!
!            = 3  error in the specification of ityp
!
!            = 4  error in the specification of nt
!
!            = 5  error in the specification of ids
!
!            = 6  error in the specification of jds
!
!            = 7  error in the specification of mdbc
!
!            = 8  error in the specification of ndbc
!
!            = 9  error in the specification of lshsgs
!
!            = 10 error in the specification of lwork
!
!
! **********************************************************************
!                                                                              
!     end of documentation for islapgs
!
! **********************************************************************
!
!
submodule(scalar_laplacian_routines) invert_scalar_laplacian_gaussian_grid_saved

contains

    module subroutine islapgs(nlat, nlon, isym, nt, xlmbda, sf, ids, jds, a, b, &
        mdab, ndab, wshsgs, lshsgs, work, lwork, pertrb, ierror)

        ! Dummy arguments
        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: nlon
        integer(ip), intent(in)  :: isym
        integer(ip), intent(in)  :: nt
        real(wp),    intent(in)  :: xlmbda(nt)
        real(wp),    intent(out) :: sf(ids, jds, nt)
        integer(ip), intent(in)  :: ids
        integer(ip), intent(in)  :: jds
        real(wp),    intent(in)  :: a(mdab, ndab, nt)
        real(wp),    intent(in)  :: b(mdab, ndab, nt)
        integer(ip), intent(in)  :: mdab
        integer(ip), intent(in)  :: ndab
        real(wp),    intent(in)  :: wshsgs(lshsgs)
        integer(ip), intent(in)  :: lshsgs
        real(wp),    intent(out) :: work(lwork)
        integer(ip), intent(in)  :: lwork
        real(wp),    intent(out) :: pertrb(nt)
        integer(ip), intent(out) :: ierror

        ! Local variables
        integer(ip) :: ia
        integer(ip) :: ib
        integer(ip) :: ierror
        integer(ip) :: ifn
        integer(ip) :: imid
        integer(ip) :: iwk
        integer(ip) :: l1
        integer(ip) :: l2
        integer(ip) :: lp
        integer(ip) :: ls
        integer(ip) :: lwk
        integer(ip) :: lwkmin
        integer(ip) :: mmax
        integer(ip) :: mn
        integer(ip) :: nln

        ! Check input arguments
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
        if ((isym == 0 .and. ids<nlat) .or. &
            (isym>0 .and. ids<imid)) return
        ierror = 6
        if (jds < nlon) return
        ierror = 7
        mmax = min(nlat, nlon/2+1)
        if (mdab < mmax) return
        ierror = 8
        if (ndab < nlat) return
        ierror = 9
        !
        !     set and verify saved work space length
        !
        imid = (nlat+1)/2
        l2 = (nlat+mod(nlat, 2))/2
        l1 = min((nlon+2)/2, nlat)
        lp=nlat*(3*(l1+l2)-2)+(l1-1)*(l2*(2*nlat-l1)-3*l1)/2+nlon+15
        if (lshsgs < lp) return
        ierror = 10
        !
        !     set and verify unsaved work space length
        !
        ls = nlat
        if (isym > 0) ls = imid
        nln = nt*ls*nlon
        mn = mmax*nlat*nt
        !     lwkmin = nln+ls*nlon+2*mn+nlat
        !     if (lwork .lt. lwkmin) return
        l2 = (nlat+1)/2
        l1 = min(nlat, nlon/2+1)
        if (isym == 0) then
            lwkmin = (nt+1)*nlat*nlon + nlat*(2*nt*l1+1)
        else
            lwkmin = (nt+1)*l2*nlon + nlat*(2*nt*l1+1)
        end if
        if (lwork < lwkmin) return
        ierror = 0

        ! Check sign of xlmbda
        if (any(xlmbda < ZERO)) then
            ierror = -1
            return
        end if

        ! Set workspace pointers
        ia = 1
        ib = ia+mn
        ifn = ib+mn
        iwk = ifn+nlat
        lwk = lwork-2*mn-nlat

        call islapgs_lower_routine(nlat, nlon, isym, nt, xlmbda, sf, ids, jds, a, b, mdab, ndab, &
            work(ia), work(ib), mmax, work(ifn), wshsgs, lshsgs, work(iwk), lwk, &
            pertrb, ierror)

    end subroutine islapgs

    subroutine islapgs_lower_routine(nlat, nlon, isym, nt, xlmbda, sf, ids, jds, a, b, &
        mdab, ndab, as, bs, mmax, fnn, wsav, lsav, wk, lwk, pertrb, ierror)

        real(wp) :: a
        real(wp) :: as
        real(wp) :: b
        real(wp) :: bs
        
        real(wp) :: fnn
        integer(ip) :: ids
        integer(ip) :: ierror
        integer(ip) :: isym
        integer(ip) :: jds
        
        integer(ip) :: lsav
        integer(ip) :: lwk
        
        integer(ip) :: mdab
        integer(ip) :: mmax
        
        integer(ip) :: ndab
        integer(ip) :: nlat
        integer(ip) :: nlon
        integer(ip) :: nt
        real(wp) :: pertrb
        real(wp) :: sf
        real(wp) :: wk
        real(wp) :: wsav
        real(wp) :: xlmbda
        dimension sf(ids, jds, nt), a(mdab, ndab, nt), b(mdab, ndab, nt)
        dimension as(mmax, nlat, nt), bs(mmax, nlat, nt), fnn(nlat)
        dimension wsav(lsav), wk(lwk), xlmbda(nt), pertrb(nt)

        call perform_setup_for_inversion(a, b, as, bs, fnn, xlmbda, pertrb)

        ! Synthesize as, bs into sf
        call shsgs(nlat, nlon, isym, nt, sf, ids, jds, as, bs, mmax, nlat, &
            wsav, lsav, wk, lwk, ierror)


    end subroutine islapgs_lower_routine

end submodule invert_scalar_laplacian_gaussian_grid_saved