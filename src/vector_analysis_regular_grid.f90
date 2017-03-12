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
submodule(vector_analysis_routines) vector_analysis_regular_grid

contains

    ! Purpose:
    !
    !     subroutine vhaec(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, &
    !                      mdab, ndab, wvhaec, ierror)
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
    !            the program that calls vhaec. if ityp <= 2 then idvw
    !            must be at least nlat.  if ityp > 2 and nlat is
    !            even then idvw must be at least nlat/2. if ityp > 2
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
    !            4*nlat*l2+3*max(l1-2, 0)*(2*nlat-l1-1)+nlon+15
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
    !
    module subroutine vhaec(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, &
        mdab, ndab, wvhaec, ierror)

        ! Dummy arguments
        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: nlon
        integer(ip), intent(in)  :: ityp
        integer(ip), intent(in)  :: nt
        real(wp),    intent(in)  :: v(idvw, jdvw, nt)
        real(wp),    intent(in)  :: w(idvw, jdvw, nt)
        integer(ip), intent(in)  :: idvw
        integer(ip), intent(in)  :: jdvw
        real(wp),    intent(out) :: br(mdab, ndab, nt)
        real(wp),    intent(out) :: bi(mdab, ndab, nt)
        real(wp),    intent(out) :: cr(mdab, ndab, nt)
        real(wp),    intent(out) :: ci(mdab, ndab, nt)
        integer(ip), intent(in)  :: mdab
        integer(ip), intent(in)  :: ndab
        real(wp),    intent(in)  :: wvhaec(:)
        integer(ip), intent(out) :: ierror

        ! Local variables
        integer(ip) :: mmax, lzz1
        integer(ip) :: labc, lnl, imid
        integer(ip) :: ist, idv, lwork

        associate (lvhaec => size(wvhaec))

            imid = (nlat+1)/2
            mmax = min(nlat, (nlon+1)/2)
            lzz1 = 2*nlat*imid
            labc = 3*(max(mmax-2, 0)*(2*nlat-mmax-1))/2

            if (nlat < 3) then
                ierror = 1
            else if (nlon < 1) then
                ierror = 2
            else if (ityp < 0 .or. ityp > 8) then
                ierror = 4
            else if (nt < 0) then
                ierror = 4
            else if ((ityp <= 2 .and. idvw < nlat) .or. &
                (ityp > 2 .and. idvw < imid)) then
                ierror = 5
            else if (jdvw < nlon) then
                ierror = 6
            else if (mdab < mmax) then
                ierror = 7
            else if (ndab < nlat) then
                ierror = 8
            else if (lvhaec < 2*(lzz1+labc)+nlon+15) then
                ierror = 9
            else
                ierror = 0
            end if

            ! Check error flag
            if (ierror /= 0) return

            ! Set required workspace size
            select case (ityp)
                case (0:2)
                    lwork = nlat*(2*nt*nlon+max(6*imid, nlon))
                case default
                    lwork = imid*(2*nt*nlon+max(6*nlat, nlon))
            end select

            select case (ityp)
                case (0:2)
                    idv = nlat
                    ist = imid
                case default
                    idv = imid
                    ist = 0
            end select

            lnl = nt*idv*nlon

            block
                real (wp)   :: work(lwork)
                integer(ip) :: iw1, iw2, iw3, iw4, iw5
                integer(ip) :: lwzvin, jw1, jw2

                iw1 = ist+1
                iw2 = lnl+1
                iw3 = iw2+ist
                iw4 = iw2+lnl
                iw5 = iw4+3*imid*nlat
                lwzvin = lzz1+labc
                jw1 = lwzvin+1
                jw2 = jw1+lwzvin
                call vhaec_lower_utility_routine(nlat, nlon, ityp, nt, imid, idvw, jdvw, &
                    v, w, mdab, ndab, br, bi, cr, ci, idv, work, work(iw1:), work(iw2:), work(iw3:), &
                    work(iw4:), work(iw5:), wvhaec, wvhaec(jw1:), wvhaec(jw2:))
            end block
        end associate

    end subroutine vhaec

    ! Purpose:
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
    !            4*nlat*l2+3*max(l1-2, 0)*(2*nlat-l1-1)+nlon+15
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
    !
    module subroutine vhaeci(nlat, nlon, wvhaec, ierror)

        ! Dummy arguments
        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: nlon
        real(wp),    intent(out) :: wvhaec(:)
        integer(ip), intent(out) :: ierror

        ! Local variables
        type(SpherepackUtility) :: util
        integer(ip) :: imid
        integer(ip) :: iw1
        integer(ip) :: iw2
        integer(ip) :: labc
        integer(ip) :: lwzvin
        integer(ip) :: lzz1
        integer(ip) :: mmax, ldwork

        associate (lvhaec => size(wvhaec))

            imid = (nlat+1)/2
            lzz1 = 2*nlat*imid
            mmax = min(nlat, (nlon+1)/2)
            labc = 3*(max(mmax-2, 0)*(2*nlat-mmax-1))/2

            ! Check error flag
            if (nlat < 3) then
                ierror = 1
            else if (nlon < 1) then
                ierror = 2
            else if (lvhaec < 2*(lzz1+labc)+nlon+15) then
                ierror = 3
            else
                ierror = 0
            end if

            ! Check error flag
            if (ierror /= 0) return

            ! Set required workspace size
            ldwork = 2*nlat+2

            block
                real(wp) :: dwork(ldwork)

                call util%zvinit(nlat, nlon, wvhaec, dwork)

                ! Set workspace index pointers
                lwzvin = lzz1+labc
                iw1 = lwzvin+1
                iw2 = iw1+lwzvin

                call util%zwinit(nlat, nlon, wvhaec(iw1:), dwork)

                call util%hfft%initialize(nlon, wvhaec(iw2:))
            end block
        end associate

    end subroutine vhaeci

    subroutine vhaec_lower_utility_routine(nlat, nlon, ityp, nt, imid, idvw, jdvw, v, w, mdab, &
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

        type(SpherepackUtility) :: util

        nlp1 = nlat+1
        tsn = TWO/nlon
        fsn = FOUR/nlon
        mlat = mod(nlat, 2)
        mlon = mod(nlon, 2)
        mmax = min(nlat, (nlon+1)/2)
        imm1 = imid
        if (mlat /= 0) imm1 = imid-1
        if (ityp > 2) goto 3
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
        goto 2
        3 do k=1, nt
            do i=1, imm1
                do j=1, nlon
                    ve(i, j, k) = fsn*v(i, j, k)
                    vo(i, j, k) = fsn*v(i, j, k)
                    we(i, j, k) = fsn*w(i, j, k)
                    wo(i, j, k) = fsn*w(i, j, k)
                end do
            end do
        end do
2       if (mlat == 0) goto 7
        do k=1, nt
            do j=1, nlon
                ve(imid, j, k) = tsn*v(imid, j, k)
                we(imid, j, k) = tsn*w(imid, j, k)
            end do
        end do
        7 do k=1, nt
            call util%hfft%forward(idv, nlon, ve(1, 1, k), idv, wrfft)
            call util%hfft%forward(idv, nlon, we(1, 1, k), idv, wrfft)
        end do
        ndo1 = nlat
        ndo2 = nlat
        if (mlat /= 0) ndo1 = nlat-1
        if (mlat == 0) ndo2 = nlat-1
        if (ityp==2 .or. ityp==5 .or. ityp==8) goto 11
        do k=1, nt
            do mp1=1, mmax
                do np1=mp1, nlat
                    br(mp1, np1, k) = ZERO
                    bi(mp1, np1, k) = ZERO
                end do
            end do
        end do
11      if (ityp==1 .or. ityp==4 .or. ityp==7) goto 13
        do k=1, nt
            do mp1=1, mmax
                do np1=mp1, nlat
                    cr(mp1, np1, k) = ZERO
                    ci(mp1, np1, k) = ZERO
                end do
            end do
        end do
13      select case (ityp)
            case (0)
                goto 1
            case (1)
                goto 100
            case (2)
                goto 200
            case (3)
                goto 300
            case (4)
                goto 400
            case (5)
                goto 500
            case (6)
                goto 600
            case (7)
                goto 700
            case (8)
                goto 800
        end select
        !
        ! case ityp=0 ,  no symmetries
        !
1       call util%zvin(0, nlat, nlon, 0, zv, iv, wzvin)
        !
        ! case m=0
        !
        do k=1, nt
            do i=1, imid
                do np1=2, ndo2, 2
                    br(1, np1, k) = br(1, np1, k)+zv(i, np1, iv)*ve(i, 1, k)
                    cr(1, np1, k) = cr(1, np1, k)-zv(i, np1, iv)*we(i, 1, k)
                end do
            end do
        end do
        do k=1, nt
            do i=1, imm1
                do np1=3, ndo1, 2
                    br(1, np1, k) = br(1, np1, k)+zv(i, np1, iv)*vo(i, 1, k)
                    cr(1, np1, k) = cr(1, np1, k)-zv(i, np1, iv)*wo(i, 1, k)
                end do
            end do
        end do
        !
        ! case m = 1 through nlat-1
        !
        if (mmax < 2) return
        do mp1=2, mmax
            m = mp1-1
            mp2 = mp1+1
            call util%zvin(0, nlat, nlon, m, zv, iv, wzvin)
            call util%zwin(0, nlat, nlon, m, zw, iw, wzwin)
            if (mp1 > ndo1) goto 17

            do k=1, nt
                do i=1, imm1
                    do np1=mp1, ndo1, 2
                        br(mp1, np1, k) = br(mp1, np1, k)+zv(i, np1, iv)*vo(i, 2*mp1-2, k) &
                            +zw(i, np1, iw)*we(i, 2*mp1-1, k)
                        bi(mp1, np1, k) = bi(mp1, np1, k)+zv(i, np1, iv)*vo(i, 2*mp1-1, k) &
                            -zw(i, np1, iw)*we(i, 2*mp1-2, k)
                        cr(mp1, np1, k) = cr(mp1, np1, k)-zv(i, np1, iv)*wo(i, 2*mp1-2, k) &
                            +zw(i, np1, iw)*ve(i, 2*mp1-1, k)
                        ci(mp1, np1, k) = ci(mp1, np1, k)-zv(i, np1, iv)*wo(i, 2*mp1-1, k) &
                            -zw(i, np1, iw)*ve(i, 2*mp1-2, k)
                    end do
                end do
            end do

            if (mlat == 0) goto 17

            do k=1, nt
                do np1=mp1, ndo1, 2
                    br(mp1, np1, k) = br(mp1, np1, k)+zw(imid, np1, iw)*we(imid, 2*mp1-1, k)
                    bi(mp1, np1, k) = bi(mp1, np1, k)-zw(imid, np1, iw)*we(imid, 2*mp1-2, k)
                    cr(mp1, np1, k) = cr(mp1, np1, k)+zw(imid, np1, iw)*ve(imid, 2*mp1-1, k)
                    ci(mp1, np1, k) = ci(mp1, np1, k)-zw(imid, np1, iw)*ve(imid, 2*mp1-2, k)
                end do
            end do

17          if (mp2 > ndo2) goto 20

            do k=1, nt
                do i=1, imm1
                    do np1=mp2, ndo2, 2
                        br(mp1, np1, k) = br(mp1, np1, k)+zv(i, np1, iv)*ve(i, 2*mp1-2, k) &
                            +zw(i, np1, iw)*wo(i, 2*mp1-1, k)
                        bi(mp1, np1, k) = bi(mp1, np1, k)+zv(i, np1, iv)*ve(i, 2*mp1-1, k) &
                            -zw(i, np1, iw)*wo(i, 2*mp1-2, k)
                        cr(mp1, np1, k) = cr(mp1, np1, k)-zv(i, np1, iv)*we(i, 2*mp1-2, k) &
                            +zw(i, np1, iw)*vo(i, 2*mp1-1, k)
                        ci(mp1, np1, k) = ci(mp1, np1, k)-zv(i, np1, iv)*we(i, 2*mp1-1, k) &
                            -zw(i, np1, iw)*vo(i, 2*mp1-2, k)
                    end do
                end do
            end do

            if (mlat == 0) goto 20

            do k=1, nt
                do np1=mp2, ndo2, 2
                    br(mp1, np1, k) = br(mp1, np1, k)+zv(imid, np1, iv)*ve(imid, 2*mp1-2, k)
                    bi(mp1, np1, k) = bi(mp1, np1, k)+zv(imid, np1, iv)*ve(imid, 2*mp1-1, k)
                    cr(mp1, np1, k) = cr(mp1, np1, k)-zv(imid, np1, iv)*we(imid, 2*mp1-2, k)
                    ci(mp1, np1, k) = ci(mp1, np1, k)-zv(imid, np1, iv)*we(imid, 2*mp1-1, k)
                end do
            end do

20      continue
        end do
        return
        !
        ! case ityp=1 ,  no symmetries but cr and ci equal zero
        !
100     call util%zvin(0, nlat, nlon, 0, zv, iv, wzvin)
        !
        ! case m=0
        !
        do k=1, nt
            do i=1, imid
                do np1=2, ndo2, 2
                    br(1, np1, k) = br(1, np1, k)+zv(i, np1, iv)*ve(i, 1, k)
                end do
            end do
        end do

        do k=1, nt
            do i=1, imm1
                do np1=3, ndo1, 2
                    br(1, np1, k) = br(1, np1, k)+zv(i, np1, iv)*vo(i, 1, k)
                end do
            end do
        end do
        !
        ! case m = 1 through nlat-1
        !
        if (mmax < 2) return
        do mp1=2, mmax
            m = mp1-1
            mp2 = mp1+1
            call util%zvin(0, nlat, nlon, m, zv, iv, wzvin)
            call util%zwin(0, nlat, nlon, m, zw, iw, wzwin)
            if (mp1 > ndo1) goto 117

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

            if (mlat == 0) goto 117

            do k=1, nt
                do np1=mp1, ndo1, 2
                    br(mp1, np1, k) = br(mp1, np1, k)+zw(imid, np1, iw)*we(imid, 2*mp1-1, k)
                    bi(mp1, np1, k) = bi(mp1, np1, k)-zw(imid, np1, iw)*we(imid, 2*mp1-2, k)
                end do
            end do

117         if (mp2 > ndo2) goto 120
            do k=1, nt
                do i=1, imm1
                    do np1=mp2, ndo2, 2
                        br(mp1, np1, k) = br(mp1, np1, k)+zv(i, np1, iv)*ve(i, 2*mp1-2, k) &
                            +zw(i, np1, iw)*wo(i, 2*mp1-1, k)
                        bi(mp1, np1, k) = bi(mp1, np1, k)+zv(i, np1, iv)*ve(i, 2*mp1-1, k) &
                            -zw(i, np1, iw)*wo(i, 2*mp1-2, k)
                    end do
                end do
            end do

            if (mlat == 0) goto 120
            do k=1, nt
                do np1=mp2, ndo2, 2
                    br(mp1, np1, k) = br(mp1, np1, k)+zv(imid, np1, iv)*ve(imid, 2*mp1-2, k)
                    bi(mp1, np1, k) = bi(mp1, np1, k)+zv(imid, np1, iv)*ve(imid, 2*mp1-1, k)
                end do
            end do
120     continue
        end do
        return
        !
        ! case ityp=2 ,  no symmetries but br and bi equal zero
        !
200     call util%zvin(0, nlat, nlon, 0, zv, iv, wzvin)
        !
        ! case m=0
        !
        do k=1, nt
            do i=1, imid
                do np1=2, ndo2, 2
                    cr(1, np1, k) = cr(1, np1, k)-zv(i, np1, iv)*we(i, 1, k)
                end do
            end do
        end do

        do k=1, nt
            do i=1, imm1
                do np1=3, ndo1, 2
                    cr(1, np1, k) = cr(1, np1, k)-zv(i, np1, iv)*wo(i, 1, k)
                end do
            end do
        end do
        !
        ! case m = 1 through nlat-1
        !
        if (mmax < 2) return
        do mp1=2, mmax
            m = mp1-1
            mp2 = mp1+1
            call util%zvin(0, nlat, nlon, m, zv, iv, wzvin)
            call util%zwin(0, nlat, nlon, m, zw, iw, wzwin)
            if (mp1 > ndo1) goto 217
            do k=1, nt
                do i=1, imm1
                    do np1=mp1, ndo1, 2
                        cr(mp1, np1, k) = cr(mp1, np1, k)-zv(i, np1, iv)*wo(i, 2*mp1-2, k) &
                            +zw(i, np1, iw)*ve(i, 2*mp1-1, k)
                        ci(mp1, np1, k) = ci(mp1, np1, k)-zv(i, np1, iv)*wo(i, 2*mp1-1, k) &
                            -zw(i, np1, iw)*ve(i, 2*mp1-2, k)
                    end do
                end do
            end do

            if (mlat == 0) goto 217
            do k=1, nt
                do np1=mp1, ndo1, 2
                    cr(mp1, np1, k) = cr(mp1, np1, k)+zw(imid, np1, iw)*ve(imid, 2*mp1-1, k)
                    ci(mp1, np1, k) = ci(mp1, np1, k)-zw(imid, np1, iw)*ve(imid, 2*mp1-2, k)
                end do
            end do

217         if (mp2 > ndo2) goto 220
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
            if (mlat == 0) goto 220
            do k=1, nt
                do np1=mp2, ndo2, 2
                    cr(mp1, np1, k) = cr(mp1, np1, k)-zv(imid, np1, iv)*we(imid, 2*mp1-2, k)
                    ci(mp1, np1, k) = ci(mp1, np1, k)-zv(imid, np1, iv)*we(imid, 2*mp1-1, k)
                end do
            end do
220     continue
        end do
        return
        !
        ! case ityp=3 ,  v even , w odd
        !
300     call util%zvin(0, nlat, nlon, 0, zv, iv, wzvin)
        !
        ! case m=0
        !
        do k=1, nt
            do i=1, imid
                do np1=2, ndo2, 2
                    br(1, np1, k) = br(1, np1, k)+zv(i, np1, iv)*ve(i, 1, k)
                end do
            end do
        end do
        do k=1, nt
            do i=1, imm1
                do np1=3, ndo1, 2
                    cr(1, np1, k) = cr(1, np1, k)-zv(i, np1, iv)*wo(i, 1, k)
                end do
            end do
        end do
        !
        ! case m = 1 through nlat-1
        !
        if (mmax < 2) return
        do mp1=2, mmax
            m = mp1-1
            mp2 = mp1+1
            call util%zvin(0, nlat, nlon, m, zv, iv, wzvin)
            call util%zwin(0, nlat, nlon, m, zw, iw, wzwin)
            if (mp1 > ndo1) goto 317
            do k=1, nt
                do i=1, imm1
                    do np1=mp1, ndo1, 2
                        cr(mp1, np1, k) = cr(mp1, np1, k)-zv(i, np1, iv)*wo(i, 2*mp1-2, k) &
                            +zw(i, np1, iw)*ve(i, 2*mp1-1, k)
                        ci(mp1, np1, k) = ci(mp1, np1, k)-zv(i, np1, iv)*wo(i, 2*mp1-1, k) &
                            -zw(i, np1, iw)*ve(i, 2*mp1-2, k)
                    end do
                end do
            end do
            if (mlat == 0) goto 317
            do k=1, nt
                do np1=mp1, ndo1, 2
                    cr(mp1, np1, k) = cr(mp1, np1, k)+zw(imid, np1, iw)*ve(imid, 2*mp1-1, k)
                    ci(mp1, np1, k) = ci(mp1, np1, k)-zw(imid, np1, iw)*ve(imid, 2*mp1-2, k)
                end do
            end do
317         if (mp2 > ndo2) goto 320
            do k=1, nt
                do i=1, imm1
                    do np1=mp2, ndo2, 2
                        br(mp1, np1, k) = br(mp1, np1, k)+zv(i, np1, iv)*ve(i, 2*mp1-2, k) &
                            +zw(i, np1, iw)*wo(i, 2*mp1-1, k)
                        bi(mp1, np1, k) = bi(mp1, np1, k)+zv(i, np1, iv)*ve(i, 2*mp1-1, k) &
                            -zw(i, np1, iw)*wo(i, 2*mp1-2, k)
                    end do
                end do
            end do
            if (mlat == 0) goto 320
            do k=1, nt
                do np1=mp2, ndo2, 2
                    br(mp1, np1, k) = br(mp1, np1, k)+zv(imid, np1, iv)*ve(imid, 2*mp1-2, k)
                    bi(mp1, np1, k) = bi(mp1, np1, k)+zv(imid, np1, iv)*ve(imid, 2*mp1-1, k)
                end do
            end do
320     continue
        end do
        return
        !
        ! case ityp=4 ,  v even, w odd, and cr and ci equal 0.
        !
400     call util%zvin(1, nlat, nlon, 0, zv, iv, wzvin)
        !
        ! case m=0
        !
        do k=1, nt
            do i=1, imid
                do np1=2, ndo2, 2
                    br(1, np1, k) = br(1, np1, k)+zv(i, np1, iv)*ve(i, 1, k)
                end do
            end do
        end do
        !
        ! case m = 1 through nlat-1
        !
        if (mmax < 2) return
        do mp1=2, mmax
            m = mp1-1
            mp2 = mp1+1
            call util%zvin(1, nlat, nlon, m, zv, iv, wzvin)
            call util%zwin(1, nlat, nlon, m, zw, iw, wzwin)
            if (mp2 > ndo2) goto 420
            do k=1, nt
                do i=1, imm1
                    do np1=mp2, ndo2, 2
                        br(mp1, np1, k) = br(mp1, np1, k)+zv(i, np1, iv)*ve(i, 2*mp1-2, k) &
                            +zw(i, np1, iw)*wo(i, 2*mp1-1, k)
                        bi(mp1, np1, k) = bi(mp1, np1, k)+zv(i, np1, iv)*ve(i, 2*mp1-1, k) &
                            -zw(i, np1, iw)*wo(i, 2*mp1-2, k)
                    end do
                end do
            end do
            if (mlat == 0) goto 420
            do k=1, nt
                do  np1=mp2, ndo2, 2
                    br(mp1, np1, k) = br(mp1, np1, k)+zv(imid, np1, iv)*ve(imid, 2*mp1-2, k)
                    bi(mp1, np1, k) = bi(mp1, np1, k)+zv(imid, np1, iv)*ve(imid, 2*mp1-1, k)
                end do
            end do
420     continue
        end do
        return
        !
        ! case ityp=5   v even, w odd, and br and bi equal zero
        !
500     call util%zvin(2, nlat, nlon, 0, zv, iv, wzvin)
        !
        ! case m=0
        !
        do k=1, nt
            do i=1, imm1
                do np1=3, ndo1, 2
                    cr(1, np1, k) = cr(1, np1, k)-zv(i, np1, iv)*wo(i, 1, k)
                end do
            end do
        end do
        !
        ! case m = 1 through nlat-1
        !
        if (mmax < 2) return
        do mp1=2, mmax
            m = mp1-1
            mp2 = mp1+1
            call util%zvin(2, nlat, nlon, m, zv, iv, wzvin)
            call util%zwin(2, nlat, nlon, m, zw, iw, wzwin)
            if (mp1 > ndo1) goto 520
            do k=1, nt
                do i=1, imm1
                    do np1=mp1, ndo1, 2
                        cr(mp1, np1, k) = cr(mp1, np1, k)-zv(i, np1, iv)*wo(i, 2*mp1-2, k) &
                            +zw(i, np1, iw)*ve(i, 2*mp1-1, k)
                        ci(mp1, np1, k) = ci(mp1, np1, k)-zv(i, np1, iv)*wo(i, 2*mp1-1, k) &
                            -zw(i, np1, iw)*ve(i, 2*mp1-2, k)
                    end do
                end do
            end do
            if (mlat == 0) goto 520
            do k=1, nt
                do np1=mp1, ndo1, 2
                    cr(mp1, np1, k) = cr(mp1, np1, k)+zw(imid, np1, iw)*ve(imid, 2*mp1-1, k)
                    ci(mp1, np1, k) = ci(mp1, np1, k)-zw(imid, np1, iw)*ve(imid, 2*mp1-2, k)
                end do
            end do
520     continue
        end do
        return
        !
        ! case ityp=6 ,  v odd , w even
        !
600     call util%zvin(0, nlat, nlon, 0, zv, iv, wzvin)
        !
        ! case m=0
        !
        do k=1, nt
            do i=1, imid
                do np1=2, ndo2, 2
                    cr(1, np1, k) = cr(1, np1, k)-zv(i, np1, iv)*we(i, 1, k)
                end do
            end do
        end do
        do k=1, nt
            do i=1, imm1
                do np1=3, ndo1, 2
                    br(1, np1, k) = br(1, np1, k)+zv(i, np1, iv)*vo(i, 1, k)
                end do
            end do
        end do
        !
        ! case m = 1 through nlat-1
        !
        if (mmax < 2) return
        do mp1=2, mmax
            m = mp1-1
            mp2 = mp1+1
            call util%zvin(0, nlat, nlon, m, zv, iv, wzvin)
            call util%zwin(0, nlat, nlon, m, zw, iw, wzwin)
            if (mp1 > ndo1) goto 617
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
            if (mlat == 0) goto 617
            do k=1, nt
                do np1=mp1, ndo1, 2
                    br(mp1, np1, k) = br(mp1, np1, k)+zw(imid, np1, iw)*we(imid, 2*mp1-1, k)
                    bi(mp1, np1, k) = bi(mp1, np1, k)-zw(imid, np1, iw)*we(imid, 2*mp1-2, k)
                end do
            end do
617         if (mp2 > ndo2) goto 620
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
            if (mlat == 0) goto 620
            do k=1, nt
                do np1=mp2, ndo2, 2
                    cr(mp1, np1, k) = cr(mp1, np1, k)-zv(imid, np1, iv)*we(imid, 2*mp1-2, k)
                    ci(mp1, np1, k) = ci(mp1, np1, k)-zv(imid, np1, iv)*we(imid, 2*mp1-1, k)
                end do
            end do
620     continue
        end do
        return
        !
        ! case ityp=7   v odd, w even, and cr and ci equal zero
        !
700     call util%zvin(2, nlat, nlon, 0, zv, iv, wzvin)
        !
        ! case m=0
        !
        do k=1, nt
            do i=1, imm1
                do np1=3, ndo1, 2
                    br(1, np1, k) = br(1, np1, k)+zv(i, np1, iv)*vo(i, 1, k)
                end do
            end do
        end do
        !
        ! case m = 1 through nlat-1
        !
        if (mmax < 2) return
        do mp1=2, mmax
            m = mp1-1
            mp2 = mp1+1
            call util%zvin(2, nlat, nlon, m, zv, iv, wzvin)
            call util%zwin(2, nlat, nlon, m, zw, iw, wzwin)
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
720     continue
        end do
        return
        !
        ! case ityp=8   v odd, w even, and both br and bi equal zero
        !
800     call util%zvin(1, nlat, nlon, 0, zv, iv, wzvin)
        !
        ! case m=0
        !
        do k=1, nt
            do i=1, imid
                do np1=2, ndo2, 2
                    cr(1, np1, k) = cr(1, np1, k)-zv(i, np1, iv)*we(i, 1, k)
                end do
            end do
        end do
        !
        ! case m = 1 through nlat-1
        !
        if (mmax < 2) return
        do mp1=2, mmax
            m = mp1-1
            mp2 = mp1+1
            call util%zvin(1, nlat, nlon, m, zv, iv, wzvin)
            call util%zwin(1, nlat, nlon, m, zw, iw, wzwin)
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
820     continue
        end do

    end subroutine vhaec_lower_utility_routine

end submodule vector_analysis_regular_grid
