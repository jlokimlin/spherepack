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
! ... file vrtes.f
!
!     this file includes documentation and code for
!     subroutine divec          i
!
! ... files which must be loaded with vrtes.f
!
!     type_SpherepackUtility.f, type_RealPeriodicFastFourierTransform.f, vhaes.f, shses.f
!
!     subroutine vrtes(nlat, nlon, isym, nt, vt, ivrt, jvrt, cr, ci, mdc, ndc, 
!    +                 wshses, lshses, work, lwork, ierror)
!
!     given the vector spherical harmonic coefficients cr and ci, precomputed
!     by subroutine vhaes for a vector field (v, w), subroutine vrtes
!     computes the vorticity of the vector field in the scalar array
!     vt.  vt(i, j) is the vorticity at the colatitude
!
!            theta(i) = (i-1)*pi/(nlat-1)
!
!     and longitude
!
!            lambda(j) = (j-1)*2*pi/nlon
!
!     on the sphere.  i.e., 
!
!            vt(i, j) =  [-dv/dlambda + d(sint*w)/dtheta]/sint
!
!     where sint = sin(theta(i)).  w is the east longitudinal and v
!     is the colatitudinal component of the vector field from which
!     cr, ci were precomputed.  required associated legendre polynomials
!     are stored rather than recomputed  as they are in subroutine vrtec.
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
!            sphere.  i.e., in the array vt(i, j) for i=1, ..., nlat and
!            j=1, ..., nlon.
!
!      = 1
!            w is antisymmetric and v is symmetric about the equator.
!            in this case the vorticity is symmetyric about the
!            equator and is computed for the northern hemisphere
!            only.  i.e., if nlat is odd the vorticity is computed
!            in the array vt(i, j) for i=1, ..., (nlat+1)/2 and for
!            j=1, ..., nlon.  if nlat is even the vorticity is computed
!            in the array vt(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
!
!      = 2
!            w is symmetric and v is antisymmetric about the equator
!            in this case the vorticity is antisymmetric about the
!            equator and is computed for the northern hemisphere
!            only.  i.e., if nlat is odd the vorticity is computed
!            in the array vt(i, j) for i=1, ..., (nlat+1)/2 and for
!            j=1, ..., nlon.  if nlat is even the vorticity is computed
!            in the array vt(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
!
!
!      nt    nt is the number of scalar and vector fields.  some
!            computational efficiency is obtained for multiple fields.
!            in the program that calls vrtes, the arrays cr, ci, and vort
!            can be three dimensional corresponding to an indexed multiple
!            vector field.  in this case multiple scalar synthesis will
!            be performed to compute the vorticity for each field.  the
!            third index is the synthesis index which assumes the values
!            k=1, ..., nt.  for a single synthesis set nt = 1.  the
!            description of the remaining parameters is simplified by
!            assuming that nt=1 or that all the arrays are two dimensional.
!
!     ivrt   the first dimension of the array vt as it appears in
!            the program that calls vrtes. if isym = 0 then ivrt
!            must be at least nlat.  if isym = 1 or 2 and nlat is
!            even then ivrt must be at least nlat/2. if isym = 1 or 2
!            and nlat is odd then ivrt must be at least (nlat+1)/2.
!
!     jvrt   the second dimension of the array vt as it appears in
!            the program that calls vrtes. jvrt must be at least nlon.
!
!    cr, ci   two or three dimensional arrays (see input parameter nt)
!            that contain vector spherical harmonic coefficients
!            of the vector field (v, w) as computed by subroutine vhaes.
!     ***    cr and ci must be computed by vhaes prior to calling
!            vrtes.
!
!      mdc   the first dimension of the arrays cr and ci as it
!            appears in the program that calls vrtes. mdc must be at
!            least min(nlat, nlon/2) if nlon is even or at least
!            min(nlat, (nlon+1)/2) if nlon is odd.
!
!      ndc   the second dimension of the arrays cr and ci as it
!            appears in the program that calls vrtes. ndc must be at
!            least nlat.
!
!   wshses   an array which must be initialized by subroutine shsesi.
!            once initialized, 
!            wshses can be used repeatedly by vrtes as long as nlon
!            and nlat remain unchanged.  wshses must not be altered
!            between calls of vrtes
!
!   lshses   the dimension of the array wshses as it appears in the
!            program that calls vrtes. define
!
!               l1 = min(nlat, (nlon+2)/2) if nlon is even or
!               l1 = min(nlat, (nlon+1)/2) if nlon is odd
!
!            and
!
!               l2 = nlat/2        if nlat is even or
!               l2 = (nlat+1)/2    if nlat is odd
!
!            then lshses must be at least
!
!               (l1*l2*(nlat+nlat-l1+1))/2+nlon+15
!
!     work   a work array that does not have to be saved.
!
!    lwork   the dimension of the array work as it appears in the
!            program that calls vrtes. define
!
!               l1 = min(nlat, nlon/2) if nlon is even or
!               l1 = min(nlat, (nlon+1)/2) if nlon is odd
!
!            and
!
!               l2 = nlat/2        if nlat is even or
!               l2 = (nlat+1)/2    if nlat is odd.
!
!            if isym = 0 then lwork must be at least
!
!               nlat*((nt+1)*nlon+2*nt*l1+1)
!
!            if isym > 0 then lwork must be at least
!
!               (nt+1)*l2*nlon+nlat*(2*nt*l1+1)
!
!
!     **************************************************************
!
!     output parameters
!
!
!     vt     a two or three dimensional array (see input parameter nt)
!            that contains the vorticity of the vector field (v, w)
!            whose coefficients cr, ci where computed by subroutine vhaes.
!            vt(i, j) is the vorticity at the colatitude point theta(i) =
!            (i-1)*pi/(nlat-1) and longitude point lambda(j) =
!            (j-1)*2*pi/nlon. the index ranges are defined above at the
!            input parameter isym.
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
!          = 9  error in the specification of lshses
!          = 10 error in the specification of lwork
! **********************************************************************
!                                                                              
!
submodule(vorticity_routines) vorticity_regular_grid_saved

contains

    module subroutine vrtes(nlat, nlon, isym, nt, vort, ivrt, jvrt, cr, ci, mdc, ndc, &
        wshses, lshses, work, lwork, ierror)

        ! Dummy arguments
        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: nlon
        integer(ip), intent(in)  :: isym
        integer(ip), intent(in)  :: nt
        real(wp),    intent(out) :: vort(ivrt, jvrt, nt)
        integer(ip), intent(in)  :: ivrt
        integer(ip), intent(in)  :: jvrt
        real(wp),    intent(in)  :: cr(mdc, ndc, nt)
        real(wp),    intent(in)  :: ci(mdc, ndc, nt)
        integer(ip), intent(in)  :: mdc
        integer(ip), intent(in)  :: ndc
        real(wp),    intent(in)  :: wshses(lshses)
        integer(ip), intent(in)  :: lshses
        real(wp),    intent(out) :: work(lwork)
        integer(ip), intent(in)  :: lwork
        integer(ip), intent(out) :: ierror

        ! Local variables
        integer(ip) :: ia
        integer(ip) :: ib
        integer(ip) :: imid
        integer(ip) :: is
        integer(ip) :: iwk
        integer(ip) :: lpimn
        integer(ip) :: ls
        integer(ip) :: lwk
        integer(ip) :: mab
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
        if (lshses < lpimn+nlon+15) return
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
        if (lwork < nln+ls*nlon+2*mn+nlat) return
        ierror = 0
        !
        !     set work space pointers
        !
        ia = 1
        ib = ia+mn
        is = ib+mn
        iwk = is+nlat
        lwk = lwork-2*mn-nlat

        call vrtes_lower_routine(nlat, nlon, isym, nt, vort, ivrt, jvrt, cr, ci, mdc, ndc, &
            work(ia), work(ib), mab, work(is), wshses, lshses, work(iwk), lwk, &
            ierror)

    end subroutine vrtes

    subroutine vrtes_lower_routine(nlat, nlon, isym, nt, vort, ivrt, jvrt, cr, ci, mdc, ndc, &
        a, b, mab, sqnn, wsav, lwsav, wk, lwk, ierror)

        real(wp) :: a
        real(wp) :: b
        real(wp) :: ci
        real(wp) :: cr
        
        integer(ip) :: ierror
        integer(ip) :: isym
        integer(ip) :: ivrt
        integer(ip) :: jvrt
        
        integer(ip) :: lwk
        integer(ip) :: lwsav
        
        integer(ip) :: mab
        integer(ip) :: mdc
        
        
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

        call perform_setup_for_vorticity(nlon, a, b, cr, ci, sqnn)

        ! Synthesize a, b into vort
        call shses(nlat, nlon, isym, nt, vort, ivrt, jvrt, a, b, &
            mab, nlat, wsav, lwsav, wk, lwk, ierror)

    end subroutine vrtes_lower_routine

end submodule vorticity_regular_grid_saved
