!
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 1998 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                      SPHEREPACK                               *
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
! ... file igrades.f
!
!     this file includes documentation and code for
!     subroutine igrades         i
!
! ... files which must be loaded with igradec.f
!
!     type_SpherepackUtility.f, type_RealPeriodicFastFourierTransform.f, shses.f, vhaes.f
!
!     subroutine igrades(nlat, nlon, isym, nt, sf, isf, jsf, br, bi, mdb, ndb, 
!    +                   wshses, lshses, work, lwork, ierror)
!
!     let br, bi, cr, ci be the vector spherical harmonic coefficients
!     precomputed by vhaes for a vector field (v, w).  let (v', w') be
!     the irrotational component of (v, w) (i.e., (v', w') is generated
!     by assuming cr, ci are zero and synthesizing br, bi with vhses).
!     then subroutine igrades computes a scalar field sf such that
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
!     at colatitude
!
!            theta(i) = (i-1)*pi/(nlat-1)
!
!     and longitude
!
!            lambda(j) = (j-1)*2*pi/nlon
!
!     where sint = sin(theta(i)).  required associated legendre polynomials
!     are stored rather than recomputed as they are in subroutine igradec.
!
!     note:  for an irrotational vector field (v, w), subroutine igrades
!     computes a scalar field whose gradient is (v, w).  in ay case, 
!     subroutine igrades "inverts" the gradient subroutine grades.
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
!            the program that calls igrades. if isym = 0 then isf
!            must be at least nlat.  if isym = 1 or 2 and nlat is
!            even then isf must be at least nlat/2. if isym = 1 or 2
!            and nlat is odd then isf must be at least (nlat+1)/2.
!
!     jsf    the second dimension of the array sf as it appears in
!            the program that calls igrades. jsf must be at least nlon.
!
!     br, bi  two or three dimensional arrays (see input parameter nt)
!            that contain vector spherical harmonic coefficients
!            of the vector field (v, w) as computed by subroutine vhaes.
!     ***    br, bi must be computed by vhaes prior to calling igrades.
!
!     mdb    the first dimension of the arrays br and bi as it appears in
!            the program that calls igrades (and vhaes). mdb must be at
!            least min(nlat, nlon/2) if nlon is even or at least
!            min(nlat, (nlon+1)/2) if nlon is odd.
!
!     ndb    the second dimension of the arrays br and bi as it appears in
!            the program that calls igrades (and vhaes). ndb must be at
!            least nlat.
!
!
!  wshses    an array which must be initialized by subroutine igradesi
!            (or equivalently by subroutine shsesi).  once initialized, 
!            wshses can be used repeatedly by igrades as long as nlon
!            and nlat remain unchanged.  wshses must not be altered
!            between calls of igrades.
!
!
!  lshses    the dimension of the array wshses as it appears in the
!            program that calls igrades. define
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
!            then lshses must be greater than or equal to
!
!               (l1*l2*(2*nlat-l1+1))/2+nlon+15
!
!     work   a work array that does not have to be saved.
!
!     lwork  the dimension of the array work as it appears in the
!            program that calls igrades. define
!
!               l2 = nlat/2                    if nlat is even or
!               l2 = (nlat+1)/2                if nlat is odd
!               l1 = min(nlat, (nlon+2)/2) if nlon is even or
!               l1 = min(nlat, (nlon+1)/2) if nlon is odd
!
!            if isym = 0 lwork must be greater than or equal to
!
!               nlat*((nt+1)*nlon+2*nt*l1+1)
!
!            if isym > 0 lwork must be greater than or equal to
!
!               (nt+1)*l2*nlon+nlat*(2*nt*l1+1)
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
!           vhaes.  sf(i, j) is given at the gaussian colatitude theta(i)
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
!           = 9  error in the specification of lshses
!           = 10 error in the specification of lwork
!
! **********************************************************************
!
submodule(gradient_routines) invert_gradient_regular_grid_saved

contains

    subroutine igrades(nlat, nlon, isym, nt, sf, isf, jsf, br, bi, mdb, ndb, &
        wshses, lshses, work, lwork, ierror)

        ! Dummy arguments
        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: nlon
        integer(ip), intent(in)  :: isym
        integer(ip), intent(in)  :: nt
        real(wp),    intent(out) :: sf(isf, jsf, nt)
        integer(ip), intent(in)  :: isf
        integer(ip), intent(in)  :: jsf
        real(wp),    intent(in)  :: br(mdb, ndb, nt)
        real(wp),    intent(in)  :: bi(mdb, ndb, nt)
        integer(ip), intent(in)  :: mdb
        integer(ip), intent(in)  :: ndb
        real(wp),    intent(in)  :: wshses(lshses)
        integer(ip), intent(in)  :: lshses
        real(wp),    intent(out) :: work(lwork)
        integer(ip), intent(in)  :: lwork
        integer(ip), intent(out) :: ierror

        ! Local variables
        integer(ip) :: ia
        integer(ip) :: ib
        integer(ip) :: imid
        integer(ip) :: iis
        integer(ip) :: iwk
        integer(ip) :: liwk
        integer(ip) :: lpimn
        integer(ip) :: ls
        integer(ip) :: lwkmin
        integer(ip) :: mab
        integer(ip) :: mmax
        integer(ip) :: mn
        integer(ip) :: nln

        ! Check calling arguments
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
        imid = (nlat+1)/2
        lpimn = (imid*mmax*(2*nlat-mmax+1))/2
        if (lshses < lpimn+nlon+15) return
        ierror = 10
        !
        !     set minimum and verify unsaved work space
        !
        ls = nlat
        if (isym > 0) ls = imid
        nln = nt*ls*nlon
        !
        !     set first dimension for a, b (as requried by shses)
        !
        mab = min(nlat, nlon/2+1)
        mn = mab*nlat*nt
        lwkmin = nln+ls*nlon+2*mn+nlat
        if (lwork < lwkmin) return
        ierror = 0
        !
        !     set work space pointers
        !
        ia = 1
        ib = ia + mn
        iis = ib + mn
        iwk = iis + nlat
        liwk = lwork-2*mn-nlat

        call igrades_lower_utility_routine(nlat, nlon, isym, nt, sf, isf, jsf, work(ia), work(ib), mab, &
            work(iis), mdb, ndb, br, bi, wshses, lshses, work(iwk), liwk, ierror)

    end subroutine igrades

    subroutine igrades_lower_utility_routine(nlat, nlon, isym, nt, sf, isf, jsf, a, b, mab, &
        sqnn, mdb, ndb, br, bi, wshses, lshses, wk, lwk, ierror)

        real(wp) :: a
        real(wp) :: b
        real(wp) :: bi
        real(wp) :: br
        
        integer(ip) :: ierror
        integer(ip) :: isf
        integer(ip) :: isym
        integer(ip) :: jsf
        
        integer(ip) :: lshses
        integer(ip) :: lwk
        
        integer(ip) :: mab
        integer(ip) :: mdb
        
        
        integer(ip) :: ndb
        integer(ip) :: nlat
        integer(ip) :: nlon
        integer(ip) :: nt
        real(wp) :: sf
        real(wp) :: sqnn
        real(wp) :: wk
        real(wp) :: wshses
        dimension sf(isf, jsf, nt)
        dimension br(mdb, ndb, nt), bi(mdb, ndb, nt), sqnn(nlat)
        dimension a(mab, nlat, nt), b(mab, nlat, nt)
        dimension wshses(lshses), wk(lwk)

        call perform_setup_for_inversion(nlon, a, b, br, bi, sqnn)

        ! Scalar synthesize a, b into sf
        call shses(nlat, nlon, isym, nt, sf, isf, jsf, a, b, mab, nlat, &
            wshses, ierror)

    end subroutine igrades_lower_utility_routine

end  submodule invert_gradient_regular_grid_saved
