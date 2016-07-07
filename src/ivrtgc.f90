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
! ... file ivrtgc.f
!
!     this file includes documentation and code for
!     subroutine ivrtgc         i
!
! ... files which must be loaded with ivrtgc.f
!
!     type_SpherepackAux.f, type_HFFTpack.f, vhsgc.f, shagc.f, gaqd.f
!
!
!     subroutine ivrtgc(nlat, nlon, isym, nt, v, w, idvw, jdvw, a, b, mdab, ndab, 
!    +                  wvhsgc, lvhsgc, work, lwork, pertrb, ierror)
!
!     given the scalar spherical harmonic coefficients a and b, precomputed
!     by subroutine shagc for a scalar array vt, subroutine ivrtgc computes
!     a divergence free vector field (v, w) whose vorticity is vt - pertrb.
!     w is the east longitude component and v is the colatitudinal component.
!     pertrb is a constant which must be subtracted from vt for (v, w) to
!     exist (see the description of pertrb below).  usually pertrb is zero
!     or small relative to vort.  the divergence of (v, w), as computed by
!     ivrtgc, is the zero scalar field.  v(i, j) and w(i, j) are the
!     colatitudinal and east longitude velocity components at gaussian
!     colatitude theta(i) (see nlat as input parameter) and longitude
!     lambda(j) = (j-1)*2*pi/nlon.  the
!
!            vorticity(v(i, j), w(i, j))
!
!         =  [-dv/dlambda + d(sint*w)/dtheta]/sint
!
!         =  vort(i, j) - pertrb
!
!     and
!
!            divergence(v(i, j), w(i, j))
!
!         =  [d(sint*v)/dtheta + dw/dlambda]/sint
!
!         =  0.0
!
!     where sint = sin(theta(i)).  required associated legendre polynomials
!     are recomputed rather than stored as they are in subroutine ivrtgs.
!
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
!            than 3. the axisymmetric case corresponds to nlon=1.
!            the efficiency of the computation is improved when nlon
!            is a product of small prime numbers.
!
!
!     isym   this has the same value as the isym that was input to
!            subroutine shagc to compute the arrays a and b.  isym
!            determines whether (v, w) are computed on the full or half
!            sphere as follows:
!
!      = 0
!            vt is not symmetric about the equator. in this case
!            the vector field (v, w) is computed on the entire sphere.
!            i.e., in the arrays  v(i, j), w(i, j) for i=1, ..., nlat and
!            j=1, ..., nlon.
!
!      = 1
!
!            vt is symmetric about the equator. in this case w is
!            antiymmetric and v is symmetric about the equator. v
!            and w are computed on the northern hemisphere only.  i.e., 
!            if nlat is odd they are computed for i=1, ..., (nlat+1)/2
!            and j=1, ..., nlon.  if nlat is even they are computed for
!            i=1, ..., nlat/2 and j=1, ..., nlon.
!
!       = 2
!
!            vt is antisymmetric about the equator. in this case w is
!            symmetric and v is antisymmetric about the equator. w
!            and v are computed on the northern hemisphere only.  i.e., 
!            if nlat is odd they are computed for i=1, ..., (nlat+1)/2
!            and j=1, ..., nlon.  if nlat is even they are computed for
!            i=1, ..., nlat/2 and j=1, ..., nlon.
!
!
!     nt     in the program that calls ivrtgc, nt is the number of vorticity
!            and vector fields.  some computational efficiency is obtained
!            for multiple fields.  the arrays a, b, v, and w can be three
!            dimensional and pertrb can be one dimensional corresponding
!            to an indexed multiple array vort.  in this case, multiple vector
!            synthesis will be performed to compute each vector field.  the
!            third index for a, b, v, w and first for pertrb is the synthesis
!            index which assumes the values k=1, ..., nt.  for a single
!            synthesis set nt=1. the description of the remaining parameters
!            is simplified by assuming that nt=1 or that a, b, v, w are two
!            dimensional and pertrb is a constant.
!
!     idvw   the first dimension of the arrays v, w as it appears in
!            the program that calls ivrtgc. if isym = 0 then idvw
!            must be at least nlat.  if isym = 1 or 2 and nlat is
!            even then idvw must be at least nlat/2. if isym = 1 or 2
!            and nlat is odd then idvw must be at least (nlat+1)/2.
!
!     jdvw   the second dimension of the arrays v, w as it appears in
!            the program that calls ivrtgc. jdvw must be at least nlon.
!
!     a, b    two or three dimensional arrays (see input parameter nt)
!            that contain scalar spherical harmonic coefficients
!            of the vorticity array vt as computed by subroutine shagc.
!     ***    a, b must be computed by shagc prior to calling ivrtgc.
!
!     mdab   the first dimension of the arrays a and b as it appears in
!            the program that calls ivrtgcs (and shagc). mdab must be at
!            least min(nlat, (nlon+2)/2) if nlon is even or at least
!            min(nlat, (nlon+1)/2) if nlon is odd.
!
!     ndab   the second dimension of the arrays a and b as it appears in
!            the program that calls ivrtgc (and shagc). ndab must be at
!            least nlat.
!
!
!  wvhsgc    an array which must be initialized by subroutine vhsgci.
!            once initialized
!            wvhsgc can be used repeatedly by ivrtgc as long as nlon
!            and nlat remain unchanged.  wvhsgs must not be altered
!            between calls of ivrtgc.
!
!
!  lvhsgc    the dimension of the array wvhsgc as it appears in the
!            program that calls ivrtgc. define
!
!               l1 = min(nlat, nlon/2) if nlon is even or
!               l1 = min(nlat, (nlon+1)/2) if nlon is odd
!
!            and
!
!               l2 = nlat/2        if nlat is even or
!               l2 = (nlat+1)/2    if nlat is odd
!
!            then lvhsgc be at least
!
!               4*nlat*l2+3*max(l1-2, 0)*(2*nlat-l1-1)+nlon+15
!
!
!     work   a work array that does not have to be saved.
!
!     lwork  the dimension of the array work as it appears in the
!            program that calls ivrtgc. define
!
!               l2 = nlat/2                    if nlat is even or
!               l2 = (nlat+1)/2                if nlat is odd
!               l1 = min(nlat, nlon/2)         if nlon is even or
!               l1 = min(nlat, (nlon+1)/2)     if nlon is odd
!
!            if isym = 0 then lwork must be at least
!
!               nlat*(2*nt*nlon+max(6*l2, nlon) + 2*nt*l1 + 1)
!
!            if isym = 1 or 2 then lwork must be at least
!
!               l2*(2*nt*nlon+max(6*nlat, nlon)) + nlat*(2*nt*l1+1)
!
!     **************************************************************
!
!     output parameters
!
!
!     v, w   two or three dimensional arrays (see input parameter nt) that
!           contain a divergence free vector field whose vorticity is
!           vt - pertrb at the gaussian colatitude point theta(i)
!           and longitude point lambda(j)=(j-1)*2*pi/nlon.  w is the east
!           longitude component and v is the colatitudinal component.  the
!           indices for v and w are defined at the input parameter isym.
!           the divergence of (v, w) is the zero scalar field.
!
!   pertrb  a nt dimensional array (see input parameter nt and assume nt=1
!           for the description that follows).  vt - pertrb is a scalar
!           field which can be the vorticity of a vector field (v, w).
!           pertrb is related to the scalar harmonic coefficients a, b
!           of vt (computed by shagc) by the formula
!
!                pertrb = a(1, 1)/(2.*sqrt(2.))
!
!           an unperturbed vt can be the vorticity of a vector field
!           only if a(1, 1) is zero.  if a(1, 1) is nonzero (flagged by
!           pertrb nonzero) then subtracting pertrb from vt yields a
!           scalar field for which a(1, 1) is zero.
!
!    ierror = 0  no errors
!           = 1  error in the specification of nlat
!           = 2  error in the specification of nlon
!           = 3  error in the specification of isym
!           = 4  error in the specification of nt
!           = 5  error in the specification of idvw
!           = 6  error in the specification of jdvw
!           = 7  error in the specification of mdab
!           = 8  error in the specification of ndab
!           = 9  error in the specification of lvhsgc
!           = 10 error in the specification of lwork
! **********************************************************************
!                                                                              
!
module module_ivrtgc

    use spherepack_precision, only: &
        wp, & ! working precision
        ip ! integer precision

    use module_vhsgc, only: &
        vhsgc

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: ivrtgc

contains

    subroutine ivrtgc(nlat, nlon, isym, nt, v, w, idvw, jdvw, a, b, mdab, ndab, &
        wvhsgc, lvhsgc, work, lwork, pertrb, ierror)

        real (wp) :: a
        real (wp) :: b
        integer (ip) :: ici
        integer (ip) :: icr
        integer (ip) :: idvw
        integer (ip) :: ierror
        integer (ip) :: imid
        integer (ip) :: is
        integer (ip) :: isym
        integer (ip) :: iwk
        integer (ip) :: jdvw
        integer (ip) :: l1
        integer (ip) :: l2
        integer (ip) :: labc
        integer (ip) :: liwk
        integer (ip) :: lvhsgc
        integer (ip) :: lwmin
        integer (ip) :: lwork
        integer (ip) :: lzz1
        integer (ip) :: mdab
        integer (ip) :: mmax
        integer (ip) :: mn
        integer (ip) :: ndab
        integer (ip) :: nlat
        integer (ip) :: nlon
        integer (ip) :: nt
        real (wp) :: pertrb
        real (wp) :: v
        real (wp) :: w
        real (wp) :: work
        real (wp) :: wvhsgc
        dimension v(idvw, jdvw, nt), w(idvw, jdvw, nt), pertrb(nt)
        dimension a(mdab, ndab, nt), b(mdab, ndab, nt)
        dimension wvhsgc(lvhsgc), work(lwork)
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
        lzz1 = 2*nlat*imid
        labc = 3*(max(mmax-2, 0)*(nlat+nlat-mmax-1))/2
        !
        !     check save work space length
        !
        l1 = min(nlat, (nlon+1)/2)
        l2 = (nlat+1)/2
        lwmin =   4*nlat*l2+3*max(l1-2, 0)*(2*nlat-l1-1)+nlon+15
        if (lvhsgc < lwmin) return
        !     if (lvhsgc .lt. 2*(lzz1+labc)+nlon+15) return
        ierror = 10
        !
        !     verify unsaved work space length
        !
        mn = mmax*nlat*nt
        if (isym /= 0  .and. lwork < &
            nlat*(2*nt*nlon+max(6*imid, nlon))+2*mn+nlat) return
        if (isym == 0  .and. lwork < &
            imid*(2*nt*nlon+max(6*nlat, nlon))+2*mn+nlat) return
        ierror = 0
        !
        !     set work space pointers
        !
        icr = 1
        ici = icr + mn
        is = ici + mn
        iwk = is + nlat
        liwk = lwork-2*mn-nlat

        call ivtgc1(nlat, nlon, isym, nt, v, w, idvw, jdvw, work(icr), work(ici), &
            mmax, work(is), mdab, ndab, a, b, wvhsgc, lvhsgc, work(iwk), &
            liwk, pertrb, ierror)

    contains

        subroutine ivtgc1(nlat, nlon, isym, nt, v, w, idvw, jdvw, cr, ci, mmax, &
            sqnn, mdab, ndab, a, b, wsav, lsav, wk, lwk, pertrb, ierror)

            real (wp) :: a
            real (wp) :: b
            real (wp) :: bi(mmax, nlat, nt)
            real (wp) :: br(mmax, nlat, nt)
            real (wp) :: ci
            real (wp) :: cr
            real (wp) :: fn
            integer (ip) :: idvw
            integer (ip) :: ierror
            integer (ip) :: isym
            integer (ip) :: ityp
            integer (ip) :: jdvw
            integer (ip) :: k
            integer (ip) :: lsav
            integer (ip) :: lwk
            integer (ip) :: m
            integer (ip) :: mdab
            integer (ip) :: mmax
            integer (ip) :: n
            integer (ip) :: ndab
            integer (ip) :: nlat
            integer (ip) :: nlon
            integer (ip) :: nt
            real (wp) :: pertrb
            real (wp) :: sqnn
            real (wp) :: v
            real (wp) :: w
            real (wp) :: wk
            real (wp) :: wsav
            dimension v(idvw, jdvw, nt), w(idvw, jdvw, nt), pertrb(nt)
            dimension cr(mmax, nlat, nt), ci(mmax, nlat, nt), sqnn(nlat)
            dimension a(mdab, ndab, nt), b(mdab, ndab, nt)
            dimension wsav(lsav), wk(lwk)
            !
            !     preset coefficient multiplyers in vector
            !
            do n=2, nlat
                fn = real(n - 1)
                sqnn(n) = sqrt(fn * (fn + 1.0))
            end do
            !
            !     compute multiple vector fields coefficients
            !
            do k=1, nt
                !
                !     set vorticity field perturbation adjustment
                !
                pertrb(k) = a(1, 1, k)/(2.*sqrt(2.))
                !
                !     preset br, bi to 0.0
                !
                do n=1, nlat
                    do m=1, mmax
                        cr(m, n, k) = 0.0
                        ci(m, n, k) = 0.0
                    end do
                end do
                !
                !     compute m=0 coefficients
                !
                do n=2, nlat
                    cr(1, n, k) = a(1, n, k)/sqnn(n)
                    ci(1, n, k) = b(1, n, k)/sqnn(n)
                end do
                !
                !     compute m>0 coefficients
                !
                do m=2, mmax
                    do n=m, nlat
                        cr(m, n, k) = a(m, n, k)/sqnn(n)
                        ci(m, n, k) = b(m, n, k)/sqnn(n)
                    end do
                end do
            end do
            !
            !     set ityp for vector synthesis with divergence=0
            !
            select case (isym)
                case (0)
                    ityp = 2
                case (1)
                    ityp = 5
                case (2)
                    ityp = 8
            end select
            !
            !     vector sythesize cr, ci into divergence free vector field (v, w)
            !
            call vhsgc(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, &
                mmax, nlat, wsav, lsav, wk, lwk, ierror)

        end subroutine ivtgc1

    end subroutine ivrtgc

end module module_ivrtgc
