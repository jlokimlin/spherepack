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
! ... file gradec.f
!
!     this file includes documentation and code for
!     subroutine gradec         i
!
! ... files which must be loaded with gradec.f
!
!     type_SpherepackUtility.f, type_RealPeriodicFastFourierTransform.f, shaec.f, vhsec.f
!
!     subroutine gradec(nlat, nlon, isym, nt, v, w, idvw, jdvw, a, b, mdab, ndab, 
!    +                  wvhsec, lvhsec, work, lwork, ierror)
!
!     given the scalar spherical harmonic coefficients a and b, precomputed
!     by subroutine shaec for a scalar field sf, subroutine gradec computes
!     an irrotational vector field (v, w) such that
!
!           gradient(sf) = (v, w).
!
!     v is the colatitudinal and w is the east longitudinal component
!     of the gradient.  i.e., 
!
!            v(i, j) = d(sf(i, j))/dtheta
!
!     and
!
!            w(i, j) = 1/sint*d(sf(i, j))/dlambda
!
!     at colatitude
!
!            theta(i) = (i-1)*pi/(nlat-1)
!
!     and longitude
!
!            lambda(j) = (j-1)*2*pi/nlon.
!
!     where sint = sin(theta(i)).  required associated legendre polynomials
!     are recomputed rather than stored as they are in subroutine grades. this
!     saves storage (compare wvhsec here and wvhses in grades) but increases
!     computational requirements.
!
!
!     input parameters
!
!
!     nlat   the number of colatitudes on the full sphere including the
!            poles. for example, nlat = 37 for a five degree grid.
!            nlat determines the grid increment in colatitude as
!            pi/(nlat-1).  if nlat is odd the equator is located at
!            grid point i=(nlat+1)/2. if nlat is even the equator is
!            located half way between points i=nlat/2 and i=nlat/2+1.
!            nlat must be at least 3. note:on the half sphere, the
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
!     isym   this has the same value as the isym that was input to
!            subroutine shaec to compute the arrays a and b from the
!            scalar field sf.  isym determines whether (v, w) are
!            computed on the full or half sphere as follows:
!
!      = 0
!
!           sf is not symmetric about the equator. in this case
!           the vector field (v, w) is computed on the entire sphere.
!           i.e., in the arrays  v(i, j), w(i, j) for i=1, ..., nlat and
!           j=1, ..., nlon.
!
!      = 1
!
!           sf is antisymmetric about the equator. in this case w is
!           antisymmetric and v is symmetric about the equator. w
!           and v are computed on the northern hemisphere only.  i.e., 
!           if nlat is odd they are computed for i=1, ..., (nlat+1)/2
!           and j=1, ..., nlon.  if nlat is even they are computed for
!           i=1, ..., nlat/2 and j=1, ..., nlon.
!
!      = 2
!
!           sf is symmetric about the equator. in this case w is
!           symmetric and v is antisymmetric about the equator. w
!           and v are computed on the northern hemisphere only.  i.e., 
!           if nlat is odd they are computed for i=1, ..., (nlat+1)/2
!           and j=1, ..., nlon.  if nlat is even they are computed for
!           i=1, ..., nlat/2 and j=1, ..., nlon.
!
!
!     nt     nt is the number of scalar and vector fields.  some
!            computational efficiency is obtained for multiple fields.
!            the arrays a, b, v, and w can be three dimensional corresponding
!            to an indexed multiple array sf.  in this case, multiple
!            vector synthesis will be performed to compute each vector
!            field.  the third index for a, b, v, and w is the synthesis
!            index which assumes the values k = 1, ..., nt.  for a single
!            synthesis set nt = 1.  the description of the remaining
!            parameters is simplified by assuming that nt=1 or that a, b, v, 
!            and w are two dimensional arrays.
!
!     idvw   the first dimension of the arrays v, w as it appears in
!            the program that calls gradec. if isym = 0 then idvw
!            must be at least nlat.  if isym = 1 or 2 and nlat is
!            even then idvw must be at least nlat/2. if isym = 1 or 2
!            and nlat is odd then idvw must be at least (nlat+1)/2.
!
!     jdvw   the second dimension of the arrays v, w as it appears in
!            the program that calls gradec. jdvw must be at least nlon.
!
!     a, b    two or three dimensional arrays (see input parameter nt)
!            that contain scalar spherical harmonic coefficients
!            of the scalar field array sf as computed by subroutine shaec.
!     ***    a, b must be computed by shaec prior to calling gradec.
!
!     mdab   the first dimension of the arrays a and b as it appears in
!            the program that calls gradec (and shaec). mdab must be at
!            least min(nlat, (nlon+2)/2) if nlon is even or at least
!            min(nlat, (nlon+1)/2) if nlon is odd.
!
!     ndab   the second dimension of the arrays a and b as it appears in
!            the program that calls gradec (and shaec). ndab must be at
!            least nlat.
!
!
!     wvhsec an array which must be initialized by subroutine vhseci.
!            once initialized, 
!            wvhsec can be used repeatedly by gradec as long as nlon
!            and nlat remain unchanged.  wvhsec must not be altered
!            between calls of gradec.
!
!
!     lvhsec the dimension of the array wvhsec as it appears in the
!            program that calls gradec. define
!
!               l1 = min(nlat, nlon/2) if nlon is even or
!               l1 = min(nlat, (nlon+1)/2) if nlon is odd
!
!            and
!
!               l2 = nlat/2        if nlat is even or
!               l2 = (nlat+1)/2    if nlat is odd.
!
!            then lvhsec must be greater than or equal to
!
!               4*nlat*l2+3*max(l1-2, 0)*(2*nlat-l1-1)+nlon+15
!
!
!     work   a work array that does not have to be saved.
!
!     lwork  the dimension of the array work as it appears in the
!            program that calls gradec. define
!
!               l1 = min(nlat, nlon/2) if nlon is even or
!               l1 = min(nlat, (nlon+1)/2) if nlon is odd
!
!            and
!
!               l2 = nlat/2                  if nlat is even or
!               l2 = (nlat+1)/2              if nlat is odd
!
!
!            if isym = 0 then lwork must be at least
!
!                nlat*(2*nt*nlon+max(6*l2, nlon)) + nlat*(2*l1*nt+1)
!
!            if isym = 1 or 2 then lwork must be at least
!
!                l2*(2*nt*nlon+max(6*nlat, nlon)) + nlat*(2*l1*nt+1)
!
!
!     **************************************************************
!
!     output parameters
!
!
!     v, w   two or three dimensional arrays (see input parameter nt) that
!           contain an irrotational vector field such that the gradient of
!           the scalar field sf is (v, w).  w(i, j) is the east longitude
!           component and v(i, j) is the colatitudinal component of velocity
!           at colatitude theta(i) = (i-1)*pi/(nlat-1) and longitude
!           lambda(j) = (j-1)*2*pi/nlon. the indices for v and w are defined
!           at the input parameter isym.  the vorticity of (v, w) is zero.
!           note that any nonzero vector field on the sphere will be
!           multiple valued at the poles [reference swarztrauber].
!
!
!  ierror   = 0  no errors
!           = 1  error in the specification of nlat
!           = 2  error in the specification of nlon
!           = 3  error in the specification of isym
!           = 4  error in the specification of nt
!           = 5  error in the specification of idvw
!           = 6  error in the specification of jdvw
!           = 7  error in the specification of mdab
!           = 8  error in the specification of ndab
!           = 9  error in the specification of lvhsec
!           = 10 error in the specification of lwork
!
!
submodule(gradient_routines) gradient_regular_grid

contains

    module subroutine gradec(nlat, nlon, isym, nt, v, w, idvw, jdvw, a, b, mdab, ndab, &
        wvhsec, lvhsec, work, lwork, ierror)

        ! Dummy arguments
        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: nlon
        integer(ip), intent(in)  :: isym
        integer(ip), intent(in)  :: nt
        real(wp),    intent(out) :: v(idvw, jdvw, nt)
        real(wp),    intent(out) :: w(idvw, jdvw, nt)
        integer(ip), intent(in)  :: idvw
        integer(ip), intent(in)  :: jdvw
        real(wp),    intent(in)  :: a(mdab, ndab, nt)
        real(wp),    intent(in)  :: b(mdab, ndab, nt)
        integer(ip), intent(in)  :: mdab
        integer(ip), intent(in)  :: ndab
        real(wp),    intent(in)  :: wvhsec(lvhsec)
        integer(ip), intent(in)  :: lvhsec
        real(wp),    intent(out) :: work(lwork)
        integer(ip), intent(in)  :: lwork
        integer(ip), intent(out) :: ierror

        ! Local variables
        integer(ip) :: ibi
        integer(ip) :: ibr
        integer(ip) :: idz
        integer(ip) :: imid
        integer(ip) :: iis
        integer(ip) :: iwk
        integer(ip) :: l1
        integer(ip) :: l2
        integer(ip) :: liwk
        integer(ip) :: lwkmin
        integer(ip) :: lwmin
        integer(ip) :: lzimn
        integer(ip) :: mmax
        integer(ip) :: mn

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
        if ((isym == 0 .and. idvw<nlat) .or. &
            (isym /= 0 .and. idvw<imid)) return
        ierror = 6
        if (jdvw < nlon) return
        ierror = 7
        mmax = min(nlat, (nlon+1)/2)
        if (mdab < min(nlat, (nlon+2)/2)) return
        ierror = 8
        if (ndab < nlat) return
        ierror = 9
        !
        !     verify minimum saved work space length
        !
        idz = (mmax*(2*nlat-mmax+1))/2
        lzimn = idz*imid
        l1 = min(nlat, (nlon+1)/2)
        l2 = (nlat+1)/2
        lwmin = 4*nlat*l2+3*max(l1-2, 0)*(2*nlat-l1-1)+nlon+15
        if (lvhsec < lwmin) return
        ierror = 10
        !
        !     verify minimum unsaved work space length
        !
        mn = mmax*nlat*nt
        if (isym == 0) then
            lwkmin = nlat*(2*nt*nlon+max(6*l2, nlon)) + nlat*(2*l1*nt+1)
        else
            lwkmin = l2*(2*nt*nlon+max(6*nlat, nlon)) + nlat*(2*l1*nt+1)
        end if
        if (lwork < lwkmin) return

        ierror = 0
        !
        !     set work space pointers
        !
        ibr = 1
        ibi = ibr + mn
        iis = ibi + mn
        iwk = iis + nlat
        liwk = lwork-2*mn-nlat
        call gradec_lower_utility_routine(nlat, nlon, isym, nt, v, w, idvw, jdvw, work(ibr), work(ibi), &
            mmax, work(iis), mdab, ndab, a, b, wvhsec, lvhsec, work(iwk), liwk, &
            ierror)

    end subroutine gradec

    subroutine gradec_lower_utility_routine(nlat, nlon, isym, nt, v, w, idvw, jdvw, br, bi, mmax, &
        sqnn, mdab, ndab, a, b, wvhsec, lvhsec, wk, lwk, ierror)

        real(wp) :: a
        real(wp) :: b
        real(wp) :: bi
        real(wp) :: br
        real(wp) :: ci(mmax, nlat, nt)
        real(wp) :: cr(mmax, nlat, nt)
        
        integer(ip) :: idvw
        integer(ip) :: ierror
        integer(ip) :: isym
        integer(ip) :: ityp
        integer(ip) :: jdvw
        
        integer(ip) :: lvhsec
        integer(ip) :: lwk
        
        integer(ip) :: mdab
        integer(ip) :: mmax
        
        integer(ip) :: ndab
        integer(ip) :: nlat
        integer(ip) :: nlon
        integer(ip) :: nt
        real(wp) :: sqnn
        real(wp) :: v
        real(wp) :: w
        real(wp) :: wk
        real(wp) :: wvhsec
        dimension v(idvw, jdvw, nt), w(idvw, jdvw, nt)
        dimension br(mmax, nlat, nt), bi(mmax, nlat, nt), sqnn(nlat)
        dimension a(mdab, ndab, nt), b(mdab, ndab, nt)
        dimension wvhsec(lvhsec), wk(lwk)

        call perform_setup_for_gradient(isym, ityp, a, b, br, bi, sqnn)

        ! Vector synthesize br, bi into (v, w) (cr, ci are dummy variables)
        call vhsec(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, &
            mmax, nlat, wvhsec, ierror)

    end subroutine gradec_lower_utility_routine

end submodule gradient_regular_grid
