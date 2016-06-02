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
! ... file vhses.f90
!
!     this file contains code and documentation for subroutines
!     vhses and vhsesi
!
! ... files which must be loaded with vhses.f90
!
!     type_SpherepackAux.f90, type_HFFTpack.f90
!
!   
!     subroutine vhses(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, 
!    +                 mdab, ndab, wvhses, lvhses, work, lwork, ierror)
!
!     subroutine vhses performs the vector spherical harmonic synthesis
!     of the arrays br, bi, cr, and ci and stores the result in the
!     arrays v and w. v(i, j) and w(i, j) are the colatitudinal 
!     (measured from the north pole) and east longitudinal components
!     respectively, located at colatitude theta(i) = (i-1)*pi/(nlat-1)
!     and longitude phi(j) = (j-1)*2*pi/nlon. the spectral
!     representation of (v, w) is given below at output parameters v, w.
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
!     ityp   = 0  no symmetries exist about the equator. the synthesis
!                 is performed on the entire sphere.  i.e. on the
!                 arrays v(i, j), w(i, j) for i=1, ..., nlat and 
!                 j=1, ..., nlon.   
!
!            = 1  no symmetries exist about the equator. the synthesis
!                 is performed on the entire sphere.  i.e. on the
!                 arrays v(i, j), w(i, j) for i=1, ..., nlat and 
!                 j=1, ..., nlon. the curl of (v, w) is zero. that is, 
!                 (d/dtheta (sin(theta) w) - dv/dphi)/sin(theta) = 0. 
!                 the coefficients cr and ci are zero.
!
!            = 2  no symmetries exist about the equator. the synthesis
!                 is performed on the entire sphere.  i.e. on the
!                 arrays v(i, j), w(i, j) for i=1, ..., nlat and 
!                 j=1, ..., nlon. the divergence of (v, w) is zero. i.e., 
!                 (d/dtheta (sin(theta) v) + dw/dphi)/sin(theta) = 0. 
!                 the coefficients br and bi are zero.
!
!            = 3  v is symmetric and w is antisymmetric about the 
!                 equator. the synthesis is performed on the northern
!                 hemisphere only.  i.e., if nlat is odd the synthesis
!                 is performed on the arrays v(i, j), w(i, j) for 
!                 i=1, ..., (nlat+1)/2 and j=1, ..., nlon. if nlat is
!                 even the synthesis is performed on the the arrays
!                 v(i, j), w(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
!
!            = 4  v is symmetric and w is antisymmetric about the 
!                 equator. the synthesis is performed on the northern
!                 hemisphere only.  i.e., if nlat is odd the synthesis
!                 is performed on the arrays v(i, j), w(i, j) for 
!                 i=1, ..., (nlat+1)/2 and j=1, ..., nlon. if nlat is
!                 even the synthesis is performed on the the arrays
!                 v(i, j), w(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
!                 the curl of (v, w) is zero. that is, 
!                 (d/dtheta (sin(theta) w) - dv/dphi)/sin(theta) = 0. 
!                 the coefficients cr and ci are zero.
!
!            = 5  v is symmetric and w is antisymmetric about the 
!                 equator. the synthesis is performed on the northern
!                 hemisphere only.  i.e., if nlat is odd the synthesis
!                 is performed on the arrays v(i, j), w(i, j) for 
!                 i=1, ..., (nlat+1)/2 and j=1, ..., nlon. if nlat is
!                 even the synthesis is performed on the the arrays
!                 v(i, j), w(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
!                 the divergence of (v, w) is zero. i.e., 
!                 (d/dtheta (sin(theta) v) + dw/dphi)/sin(theta) = 0. 
!                 the coefficients br and bi are zero.
!
!            = 6  v is antisymmetric and w is symmetric about the 
!                 equator. the synthesis is performed on the northern
!                 hemisphere only.  i.e., if nlat is odd the synthesis
!                 is performed on the arrays v(i, j), w(i, j) for 
!                 i=1, ..., (nlat+1)/2 and j=1, ..., nlon. if nlat is
!                 even the synthesis is performed on the the arrays
!                 v(i, j), w(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
!
!            = 7  v is antisymmetric and w is symmetric about the 
!                 equator. the synthesis is performed on the northern
!                 hemisphere only.  i.e., if nlat is odd the synthesis
!                 is performed on the arrays v(i, j), w(i, j) for 
!                 i=1, ..., (nlat+1)/2 and j=1, ..., nlon. if nlat is
!                 even the synthesis is performed on the the arrays
!                 v(i, j), w(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
!                 the curl of (v, w) is zero. that is, 
!                 (d/dtheta (sin(theta) w) - dv/dphi)/sin(theta) = 0. 
!                 the coefficients cr and ci are zero.
!
!            = 8  v is antisymmetric and w is symmetric about the 
!                 equator. the synthesis is performed on the northern
!                 hemisphere only.  i.e., if nlat is odd the synthesis
!                 is performed on the arrays v(i, j), w(i, j) for 
!                 i=1, ..., (nlat+1)/2 and j=1, ..., nlon. if nlat is
!                 even the synthesis is performed on the the arrays
!                 v(i, j), w(i, j) for i=1, ..., nlat/2 and j=1, ..., nlon.
!                 the divergence of (v, w) is zero. i.e., 
!                 (d/dtheta (sin(theta) v) + dw/dphi)/sin(theta) = 0. 
!                 the coefficients br and bi are zero.
!
!
!     nt     the number of syntheses.  in the program that calls vhses, 
!            the arrays v, w, br, bi, cr, and ci can be three dimensional
!            in which case multiple syntheses will be performed.
!            the third index is the synthesis index which assumes the 
!            values k=1, ..., nt.  for a single synthesis set nt=1. the
!            discription of the remaining parameters is simplified
!            by assuming that nt=1 or that all the arrays are two
!            dimensional.
!
!     idvw   the first dimension of the arrays v, w as it appears in
!            the program that calls vhaes. if ityp .le. 2 then idvw
!            must be at least nlat.  if ityp .gt. 2 and nlat is
!            even then idvw must be at least nlat/2. if ityp .gt. 2
!            and nlat is odd then idvw must be at least (nlat+1)/2.
!
!     jdvw   the second dimension of the arrays v, w as it appears in
!            the program that calls vhses. jdvw must be at least nlon.
!
!     br, bi  two or three dimensional arrays (see input parameter nt)
!     cr, ci  that contain the vector spherical harmonic coefficients
!            in the spectral representation of v(i, j) and w(i, j) given
!            below at the discription of output parameters v and w.
!
!     mdab   the first dimension of the arrays br, bi, cr, and ci as it
!            appears in the program that calls vhses. mdab must be at
!            least min(nlat, nlon/2) if nlon is even or at least
!            min(nlat, (nlon+1)/2) if nlon is odd.
!
!     ndab   the second dimension of the arrays br, bi, cr, and ci as it
!            appears in the program that calls vhses. ndab must be at
!            least nlat.
!
!     wvhses an array which must be initialized by subroutine vhsesi.
!            once initialized, wvhses can be used repeatedly by vhses
!            as long as nlon and nlat remain unchanged.  wvhses must
!            not be altered between calls of vhses.
!
!     lvhses the dimension of the array wvhses as it appears in the
!            program that calls vhses. define
!
!               l1 = min(nlat, nlon/2) if nlon is even or
!               l1 = min(nlat, (nlon+1)/2) if nlon is odd
!
!            and
!
!               l2 = nlat/2        if nlat is even or
!               l2 = (nlat+1)/2    if nlat is odd
!
!            then lvhses must be at least
!
!                 l1*l2*(nlat+nlat-l1+1)+nlon+15
!
!
!     work   a work array that does not have to be saved.
!
!     lwork  the dimension of the array work as it appears in the
!            program that calls vhses. define
!
!               l2 = nlat/2        if nlat is even or
!               l2 = (nlat+1)/2    if nlat is odd
!
!            if ityp .le. 2 then lwork must be at least
!
!                       (2*nt+1)*nlat*nlon
!
!            if ityp .gt. 2 then lwork must be at least
!
!                        (2*nt+1)*l2*nlon 
!
!     **************************************************************
!
!     output parameters
!
!     v, w    two or three dimensional arrays (see input parameter nt)
!            in which the synthesis is stored. v is the colatitudinal
!            component and w is the east longitudinal component. 
!            v(i, j), w(i, j) contain the components at colatitude
!            theta(i) = (i-1)*pi/(nlat-1) and longitude phi(j) =
!            (j-1)*2*pi/nlon. the index ranges are defined above at
!            the input parameter ityp. v and w are computed from the 
!            formulas given below
!
!
!     define
!
!     1.  theta is colatitude and phi is east longitude
!
!     2.  the normalized associated legendre funnctions
!
!         pbar(m, n, theta) = sqrt((2*n+1)*factorial(n-m)
!                        /(2*factorial(n+m)))*sin(theta)**m/(2**n*
!                        factorial(n)) times the (n+m)th derivative
!                        of (x**2-1)**n with respect to x=cos(theta)
!
!     3.  vbar(m, n, theta) = the derivative of pbar(m, n, theta) with
!                           respect to theta divided by the square
!                           root of n(n+1).
!
!         vbar(m, n, theta) is more easily computed in the form
!
!         vbar(m, n, theta) = (sqrt((n+m)*(n-m+1))*pbar(m-1, n, theta)
!         -sqrt((n-m)*(n+m+1))*pbar(m+1, n, theta))/(2*sqrt(n*(n+1)))
!
!     4.  wbar(m, n, theta) = m/(sin(theta))*pbar(m, n, theta) divided
!                           by the square root of n(n+1).
!
!         wbar(m, n, theta) is more easily computed in the form
!
!         wbar(m, n, theta) = sqrt((2n+1)/(2n-1))*(sqrt((n+m)*(n+m-1))
!         *pbar(m-1, n-1, theta)+sqrt((n-m)*(n-m-1))*pbar(m+1, n-1, theta))
!         /(2*sqrt(n*(n+1)))
!
!
!    the colatitudnal dependence of the normalized surface vector
!                spherical harmonics are defined by
!
!     5.    bbar(m, n, theta) = (vbar(m, n, theta), i*wbar(m, n, theta))
!
!     6.    cbar(m, n, theta) = (i*wbar(m, n, theta), -vbar(m, n, theta))
!
!
!    the coordinate to index mappings 
!
!     7.   theta(i) = (i-1)*pi/(nlat-1) and phi(j) = (j-1)*2*pi/nlon
!
!    
!     the maximum (plus one) longitudinal wave number
!
!     8.     mmax = min(nlat, nlon/2) if nlon is even or
!            mmax = min(nlat, (nlon+1)/2) if nlon is odd.
!
!    if we further define the output vector as
!
!     9.    h(i, j) = (v(i, j), w(i, j))
!
!    and the complex coefficients
!
!     10.   b(m, n) = cmplx(br(m+1, n+1), bi(m+1, n+1))
!
!     11.   c(m, n) = cmplx(cr(m+1, n+1), ci(m+1, n+1))
!
!
!    then for i=1, ..., nlat and  j=1, ..., nlon
!
!        the expansion for real h(i, j) takes the form
!
!     h(i, j) = the sum from n=1 to n=nlat-1 of the real part of
!
!         0.5_wp*(b(0, n)*bbar(0, n, theta(i))+c(0, n)*cbar(0, n, theta(i)))
!
!     plus the sum from m=1 to m=mmax-1 of the sum from n=m to 
!     n=nlat-1 of the real part of
!
!              b(m, n)*bbar(m, n, theta(i))*exp(i*m*phi(j))
!             +c(m, n)*cbar(m, n, theta(i))*exp(i*m*phi(j))
!
!   *************************************************************
!
!   in terms of real variables this expansion takes the form
!
!             for i=1, ..., nlat and  j=1, ..., nlon
!
!     v(i, j) = the sum from n=1 to n=nlat-1 of
!
!               0.5_wp*br(1, n+1)*vbar(0, n, theta(i))
!
!     plus the sum from m=1 to m=mmax-1 of the sum from n=m to 
!     n=nlat-1 of the real part of
!
!       (br(m+1, n+1)*vbar(m, n, theta(i))-ci(m+1, n+1)*wbar(m, n, theta(i)))
!                                          *cos(m*phi(j))
!      -(bi(m+1, n+1)*vbar(m, n, theta(i))+cr(m+1, n+1)*wbar(m, n, theta(i)))
!                                          *sin(m*phi(j))
!
!    and for i=1, ..., nlat and  j=1, ..., nlon
!
!     w(i, j) = the sum from n=1 to n=nlat-1 of
!
!              -0.5_wp*cr(1, n+1)*vbar(0, n, theta(i))
!
!     plus the sum from m=1 to m=mmax-1 of the sum from n=m to
!     n=nlat-1 of the real part of
!
!      -(cr(m+1, n+1)*vbar(m, n, theta(i))+bi(m+1, n+1)*wbar(m, n, theta(i)))
!                                          *cos(m*phi(j))
!      +(ci(m+1, n+1)*vbar(m, n, theta(i))-br(m+1, n+1)*wbar(m, n, theta(i)))
!                                          *sin(m*phi(j))
!
!
!      br(m+1, nlat), bi(m+1, nlat), cr(m+1, nlat), and ci(m+1, nlat) are
!      assumed zero for m even.
!
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
!            = 9  error in the specification of lvhses
!            = 10 error in the specification of lwork
!
! ************************************************************
!
!     subroutine vhsesi(nlat, nlon, wvhses, lvhses, work, lwork, dwork, 
!    +                  ldwork, ierror)
!
!     subroutine vhsesi initializes the array wvhses which can then be
!     used repeatedly by subroutine vhses until nlat or nlon is changed.
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
!     lvhses the dimension of the array wvhses as it appears in the
!            program that calls vhses. define
!
!               l1 = min(nlat, nlon/2) if nlon is even or
!               l1 = min(nlat, (nlon+1)/2) if nlon is odd
!
!            and
!
!               l2 = nlat/2        if nlat is even or
!               l2 = (nlat+1)/2    if nlat is odd
!
!            then lvhses must be at least
!
!                  l1*l2*(nlat+nlat-l1+1)+nlon+15
!
!
!     work   a work array that does not have to be saved.
!
!     lwork  the dimension of the array work as it appears in the
!            program that calls vhses. lwork must be at least
!
!              3*(max(l1-2, 0)*(nlat+nlat-l1-1))/2+5*l2*nlat
!
!     dwork  an unsaved real work space
!
!     ldwork the length of the array dwork as it appears in the
!            program that calls vhsesi.  ldwork must be at least
!            2*(nlat+1)
!
!
!
!     output parameters
!
!     wvhses an array which is initialized for use by subroutine vhses.
!            once initialized, wvhses can be used repeatedly by vhses
!            as long as nlat or nlon remain unchanged.  wvhses must not
!            be altered between calls of vhses.
!
!
!     ierror = 0  no errors
!            = 1  error in the specification of nlat
!            = 2  error in the specification of nlon
!            = 3  error in the specification of lvhses
!            = 4  error in the specification of lwork
!            = 5  error in the specification of ldwork
!
!
module module_vhses

    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32

    use type_HFFTpack, only: &
        HFFTpack

    use type_SpherepackAux, only: &
        SpherepackAux

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    public :: vhses
    public :: vhsesi

contains

    subroutine vhses(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, &
        mdab, ndab, wvhses, lvhses, work, lwork, ierror)
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        integer (ip), intent (in)  :: nlat
        integer (ip), intent (in)  :: nlon
        integer (ip), intent (in)  :: ityp
        integer (ip), intent (in)  :: nt
        real (wp),    intent (out) :: v(idvw, jdvw,*)
        real (wp),    intent (out) :: w(idvw, jdvw,*)
        integer (ip), intent (in)  :: idvw
        integer (ip), intent (in)  :: jdvw
        real (wp),    intent (in)  :: br(mdab, ndab,*)
        real (wp),    intent (in)  :: bi(mdab, ndab,*)
        real (wp),    intent (in)  :: cr(mdab, ndab,*)
        real (wp),    intent (in)  :: ci(mdab, ndab,*)
        integer (ip), intent (in)  :: mdab
        integer (ip), intent (in)  :: ndab
        real (wp),    intent (in)  :: wvhses(lvhses)
        integer (ip), intent (in)  :: lvhses
        real (wp),    intent (out) :: work(lwork)
        integer (ip), intent (in)  :: lwork
        integer (ip), intent (out) :: ierror
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip) :: idv, idz, imid, ist, lnl, lzimn, mmax
        integer (ip) :: workspace_indices(6)
        !----------------------------------------------------------------------

        imid = (nlat+1)/2
        mmax = min(nlat, (nlon+1)/2)
        idz = (mmax*(2*nlat-mmax+1))/2
        lzimn = idz*imid

        select case (ityp)
            case (:2)
                idv = nlat
                ist = imid
            case default
                idv = imid
                ist = 0
        end select

        lnl = nt*idv*nlon

        !
        !==> Check validity of input arguments
        !
        if (nlat < 3) then
            ierror = 1
            return
        else if (nlon < 1) then
            ierror = 2
            return
        else if (ityp < 0 .or. ityp > 8) then
            ierror = 3
            return
        else if (nt < 0) then
            ierror = 4
            return
        else if ((ityp <= 2 .and. idvw < nlat) .or. (ityp > 2 .and. idvw < imid)) then
            ierror = 5
            return
        else if (jdvw < nlon) then
            ierror = 6
            return
        else if (mdab < mmax) then
            ierror = 7
            return
        else if (ndab < nlat) then
            ierror = 8
            return
        else if (lvhses < 2*lzimn+nlon+15) then
            ierror = 9
            return
        else if (lwork < 2*lnl+idv*nlon) then
            ierror = 10
            return
        else
            ierror = 0
        end if

        !
        !==> Set workspace indices
        !
        workspace_indices = get_workspace_indices(ist, lnl, lzimn)

        associate( &
            iw1 => workspace_indices(1), &
            iw2 => workspace_indices(2), &
            iw3 => workspace_indices(3), &
            iw4 => workspace_indices(4), &
            jw1 => workspace_indices(5), &
            jw2 => workspace_indices(6) &
            )

            call vhses1(nlat, nlon, ityp, nt, imid, idvw, jdvw, v, w, mdab, ndab, &
                br, bi, cr, ci, idv, work, work(iw1), work(iw2), work(iw3), &
                work(iw4), idz, wvhses, wvhses(jw1), wvhses(jw2))

        end associate


    contains


        pure function get_workspace_indices(ist, lnl, lzimn) result (return_value)
            !----------------------------------------------------------------------
            ! Dictionary: calling arguments
            !----------------------------------------------------------------------
            integer (ip), intent (in)  :: ist
            integer (ip), intent (in)  :: lnl
            integer (ip), intent (in)  :: lzimn
            integer (ip)               :: return_value(6)
            !----------------------------------------------------------------------

            associate( i => return_value )

                i(1) = ist+1
                i(2) = lnl+1
                i(3) = i(2)+ist
                i(4) = i(2)+lnl
                i(5) = lzimn+1
                i(6) = i(5)+lzimn

            end associate

        end function get_workspace_indices

        subroutine vhses1(nlat, nlon, ityp, nt, imid, idvw, jdvw, v, w, mdab, &
            ndab, br, bi, cr, ci, idv, ve, vo, we, wo, work, idz, vb, wb, wrfft)
            !----------------------------------------------------------------------
            ! Dictionary: calling arguments
            !----------------------------------------------------------------------
            integer (ip), intent (in)  :: nlat
            integer (ip), intent (in)  :: nlon
            integer (ip), intent (in)  :: ityp
            integer (ip), intent (in)  :: nt
            integer (ip), intent (in)  :: imid
            integer (ip), intent (in)  :: idvw
            integer (ip), intent (in)  :: jdvw
            real (wp),    intent (out) :: v(idvw, jdvw,*)
            real (wp),    intent (out) :: w(idvw, jdvw,*)
            integer (ip), intent (in)  :: mdab
            integer (ip), intent (in)  :: ndab
            real (wp),    intent (in)  :: br(mdab, ndab,* )
            real (wp),    intent (in)  :: bi(mdab, ndab, *)
            real (wp),    intent (in)  :: cr(mdab, ndab, *)
            real (wp),    intent (in)  :: ci(mdab, ndab, *)
            integer (ip), intent (in)  :: idv
            real (wp),    intent (out)  :: ve(idv, nlon, *)
            real (wp),    intent (out)  :: vo(idv, nlon, *)
            real (wp),    intent (out)  :: we(idv, nlon, *)
            real (wp),    intent (out)  :: wo(idv, nlon, *)
            real (wp),    intent (out)  :: work(*)
            integer (ip), intent (in)  :: idz
            real (wp),    intent (in)  :: vb(imid, *)
            real (wp),    intent (in)  :: wb(imid, *)
            real (wp),    intent (in)  :: wrfft(*)
            !----------------------------------------------------------------------
            ! Dictionary: local variables
            !----------------------------------------------------------------------
            integer (ip)    :: i, imm1, j, k, m, mb, mlat, mlon, mmax, mn
            integer (ip)    :: mp1, mp2, ndo1, ndo2, nlp1, np1
            type (HFFTpack) :: hfft
            !----------------------------------------------------------------------

            nlp1 = nlat+1
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

            do k=1, nt
                ve(:,:, k) = 0.0_wp
                we(:,:, k) = 0.0_wp
            end do

            case_block: block
                select case (ityp)
                    case (0)
                        !
                        !     case ityp=0   no symmetries
                        !
                        !     case m = 0
                        !
                        do k=1, nt
                            do np1=2, ndo2, 2
                                do i=1, imid
                                    ve(i, 1, k)=ve(i, 1, k)+br(1, np1, k)*vb(i, np1)
                                    we(i, 1, k)=we(i, 1, k)-cr(1, np1, k)*vb(i, np1)
                                end do
                            end do
                        end do
                        do k=1, nt
                            do np1=3, ndo1, 2
                                do i=1, imm1
                                    vo(i, 1, k)=vo(i, 1, k)+br(1, np1, k)*vb(i, np1)
                                    wo(i, 1, k)=wo(i, 1, k)-cr(1, np1, k)*vb(i, np1)
                                end do
                            end do
                        end do
                        !
                        !     case m = 1 through nlat-1
                        !
                        if (mmax < 2) exit case_block

                        do mp1=2, mmax
                            m = mp1-1
                            mb = m*(nlat-1)-(m*(m-1))/2
                            mp2 = mp1+1

                            if (mp1 <= ndo1) then
                                do k=1, nt
                                    do np1=mp1, ndo1, 2
                                        mn = mb+np1
                                        do i=1, imm1
                                            vo(i, 2*mp1-2, k) = vo(i, 2*mp1-2, k)+br(mp1, np1, k)*vb(i, mn)
                                            ve(i, 2*mp1-2, k) = ve(i, 2*mp1-2, k)-ci(mp1, np1, k)*wb(i, mn)
                                            vo(i, 2*mp1-1, k) = vo(i, 2*mp1-1, k)+bi(mp1, np1, k)*vb(i, mn)
                                            ve(i, 2*mp1-1, k) = ve(i, 2*mp1-1, k)+cr(mp1, np1, k)*wb(i, mn)
                                            wo(i, 2*mp1-2, k) = wo(i, 2*mp1-2, k)-cr(mp1, np1, k)*vb(i, mn)
                                            we(i, 2*mp1-2, k) = we(i, 2*mp1-2, k)-bi(mp1, np1, k)*wb(i, mn)
                                            wo(i, 2*mp1-1, k) = wo(i, 2*mp1-1, k)-ci(mp1, np1, k)*vb(i, mn)
                                            we(i, 2*mp1-1, k) = we(i, 2*mp1-1, k)+br(mp1, np1, k)*wb(i, mn)
                                        end do

                                        if (mlat /= 0) then
                                            ve(imid, 2*mp1-2, k) = ve(imid, 2*mp1-2, k) &
                                                -ci(mp1, np1, k)*wb(imid, mn)
                                            ve(imid, 2*mp1-1, k) = ve(imid, 2*mp1-1, k) &
                                                +cr(mp1, np1, k)*wb(imid, mn)
                                            we(imid, 2*mp1-2, k) = we(imid, 2*mp1-2, k) &
                                                -bi(mp1, np1, k)*wb(imid, mn)
                                            we(imid, 2*mp1-1, k) = we(imid, 2*mp1-1, k) &
                                                +br(mp1, np1, k)*wb(imid, mn)
                                        end if
                                    end do
                                end do
                            end if

                            if (mp2 <= ndo2) then
                                do k=1, nt
                                    do np1=mp2, ndo2, 2
                                        mn = mb+np1
                                        do i=1, imm1
                                            ve(i, 2*mp1-2, k) = ve(i, 2*mp1-2, k)+br(mp1, np1, k)*vb(i, mn)
                                            vo(i, 2*mp1-2, k) = vo(i, 2*mp1-2, k)-ci(mp1, np1, k)*wb(i, mn)
                                            ve(i, 2*mp1-1, k) = ve(i, 2*mp1-1, k)+bi(mp1, np1, k)*vb(i, mn)
                                            vo(i, 2*mp1-1, k) = vo(i, 2*mp1-1, k)+cr(mp1, np1, k)*wb(i, mn)
                                            we(i, 2*mp1-2, k) = we(i, 2*mp1-2, k)-cr(mp1, np1, k)*vb(i, mn)
                                            wo(i, 2*mp1-2, k) = wo(i, 2*mp1-2, k)-bi(mp1, np1, k)*wb(i, mn)
                                            we(i, 2*mp1-1, k) = we(i, 2*mp1-1, k)-ci(mp1, np1, k)*vb(i, mn)
                                            wo(i, 2*mp1-1, k) = wo(i, 2*mp1-1, k)+br(mp1, np1, k)*wb(i, mn)
                                        end do
                                        if (mlat /= 0) then
                                            ve(imid, 2*mp1-2, k) = ve(imid, 2*mp1-2, k) &
                                                +br(mp1, np1, k)*vb(imid, mn)
                                            ve(imid, 2*mp1-1, k) = ve(imid, 2*mp1-1, k) &
                                                +bi(mp1, np1, k)*vb(imid, mn)
                                            we(imid, 2*mp1-2, k) = we(imid, 2*mp1-2, k) &
                                                -cr(mp1, np1, k)*vb(imid, mn)
                                            we(imid, 2*mp1-1, k) = we(imid, 2*mp1-1, k) &
                                                -ci(mp1, np1, k)*vb(imid, mn)
                                        end if
                                    end do
                                end do
                            end if
                        end do
                    case (1)
                        !
                        !     case ityp=1   no symmetries,  cr and ci equal zero
                        !
                        !     case m = 0
                        !
                        do k=1, nt
                            do np1=2, ndo2, 2
                                do i=1, imid
                                    ve(i, 1, k)=ve(i, 1, k)+br(1, np1, k)*vb(i, np1)
                                end do
                            end do
                        end do
                        do k=1, nt
                            do np1=3, ndo1, 2
                                do i=1, imm1
                                    vo(i, 1, k)=vo(i, 1, k)+br(1, np1, k)*vb(i, np1)
                                end do
                            end do
                        end do
                        !
                        !     case m = 1 through nlat-1
                        !
                        if (mmax < 2) exit case_block

                        do mp1=2, mmax
                            m = mp1-1
                            mb = m*(nlat-1)-(m*(m-1))/2
                            mp2 = mp1+1
                            if (mp1 <= ndo1) then
                                do k=1, nt
                                    do np1=mp1, ndo1, 2
                                        mn = mb+np1
                                        do i=1, imm1
                                            vo(i, 2*mp1-2, k) = vo(i, 2*mp1-2, k)+br(mp1, np1, k)*vb(i, mn)
                                            vo(i, 2*mp1-1, k) = vo(i, 2*mp1-1, k)+bi(mp1, np1, k)*vb(i, mn)
                                            we(i, 2*mp1-2, k) = we(i, 2*mp1-2, k)-bi(mp1, np1, k)*wb(i, mn)
                                            we(i, 2*mp1-1, k) = we(i, 2*mp1-1, k)+br(mp1, np1, k)*wb(i, mn)
                                        end do
                                        if (mlat /= 0) then
                                            we(imid, 2*mp1-2, k) = we(imid, 2*mp1-2, k) &
                                                -bi(mp1, np1, k)*wb(imid, mn)
                                            we(imid, 2*mp1-1, k) = we(imid, 2*mp1-1, k) &
                                                +br(mp1, np1, k)*wb(imid, mn)
                                        end if
                                    end do
                                end do
                            end if

                            if (mp2 <= ndo2) then
                                do k=1, nt
                                    do np1=mp2, ndo2, 2
                                        mn = mb+np1
                                        do i=1, imm1
                                            ve(i, 2*mp1-2, k) = ve(i, 2*mp1-2, k)+br(mp1, np1, k)*vb(i, mn)
                                            ve(i, 2*mp1-1, k) = ve(i, 2*mp1-1, k)+bi(mp1, np1, k)*vb(i, mn)
                                            wo(i, 2*mp1-2, k) = wo(i, 2*mp1-2, k)-bi(mp1, np1, k)*wb(i, mn)
                                            wo(i, 2*mp1-1, k) = wo(i, 2*mp1-1, k)+br(mp1, np1, k)*wb(i, mn)
                                        end do
                                        if (mlat /= 0) then
                                            ve(imid, 2*mp1-2, k) = ve(imid, 2*mp1-2, k) &
                                                +br(mp1, np1, k)*vb(imid, mn)
                                            ve(imid, 2*mp1-1, k) = ve(imid, 2*mp1-1, k) &
                                                +bi(mp1, np1, k)*vb(imid, mn)
                                        end if
                                    end do
                                end do
                            end if
                        end do
                    case (2)
                        !
                        !     case ityp=2   no symmetries,  br and bi are equal to zero
                        !
                        !     case m = 0
                        !
                        do k=1, nt
                            do np1=2, ndo2, 2
                                do i=1, imid
                                    we(i, 1, k)=we(i, 1, k)-cr(1, np1, k)*vb(i, np1)
                                end do
                            end do
                        end do

                        do k=1, nt
                            do np1=3, ndo1, 2
                                do i=1, imm1
                                    wo(i, 1, k)=wo(i, 1, k)-cr(1, np1, k)*vb(i, np1)
                                end do
                            end do
                        end do
                        !
                        !     case m = 1 through nlat-1
                        !
                        if (mmax < 2) exit case_block

                        do mp1=2, mmax
                            m = mp1-1
                            mb = m*(nlat-1)-(m*(m-1))/2
                            mp2 = mp1+1
                            if (mp1 <= ndo1) then
                                do k=1, nt
                                    do np1=mp1, ndo1, 2
                                        mn = mb+np1
                                        do i=1, imm1
                                            ve(i, 2*mp1-2, k) = ve(i, 2*mp1-2, k)-ci(mp1, np1, k)*wb(i, mn)
                                            ve(i, 2*mp1-1, k) = ve(i, 2*mp1-1, k)+cr(mp1, np1, k)*wb(i, mn)
                                            wo(i, 2*mp1-2, k) = wo(i, 2*mp1-2, k)-cr(mp1, np1, k)*vb(i, mn)
                                            wo(i, 2*mp1-1, k) = wo(i, 2*mp1-1, k)-ci(mp1, np1, k)*vb(i, mn)
                                        end do

                                        if (mlat /= 0) then
                                            ve(imid, 2*mp1-2, k) = ve(imid, 2*mp1-2, k) &
                                                -ci(mp1, np1, k)*wb(imid, mn)
                                            ve(imid, 2*mp1-1, k) = ve(imid, 2*mp1-1, k) &
                                                +cr(mp1, np1, k)*wb(imid, mn)
                                        end if
                                    end do
                                end do
                            end if

                            if (mp2 <= ndo2) then
                                do k=1, nt
                                    do np1=mp2, ndo2, 2
                                        mn = mb+np1
                                        do i=1, imm1
                                            vo(i, 2*mp1-2, k) = vo(i, 2*mp1-2, k)-ci(mp1, np1, k)*wb(i, mn)
                                            vo(i, 2*mp1-1, k) = vo(i, 2*mp1-1, k)+cr(mp1, np1, k)*wb(i, mn)
                                            we(i, 2*mp1-2, k) = we(i, 2*mp1-2, k)-cr(mp1, np1, k)*vb(i, mn)
                                            we(i, 2*mp1-1, k) = we(i, 2*mp1-1, k)-ci(mp1, np1, k)*vb(i, mn)
                                        end do
                                        if (mlat /= 0) then
                                            we(imid, 2*mp1-2, k) = we(imid, 2*mp1-2, k) &
                                                -cr(mp1, np1, k)*vb(imid, mn)
                                            we(imid, 2*mp1-1, k) = we(imid, 2*mp1-1, k) &
                                                -ci(mp1, np1, k)*vb(imid, mn)
                                        end if
                                    end do
                                end do
                            end if
                        end do
                    case (3)
                         !
                         !     case ityp=3   v even,  w odd
                         !
                         !     case m = 0
                         !
                        do k=1, nt
                            do np1=2, ndo2, 2
                                do i=1, imid
                                    ve(i, 1, k)=ve(i, 1, k)+br(1, np1, k)*vb(i, np1)
                                end do
                            end do
                        end do

                        do k=1, nt
                            do np1=3, ndo1, 2
                                do i=1, imm1
                                    wo(i, 1, k)=wo(i, 1, k)-cr(1, np1, k)*vb(i, np1)
                                end do
                            end do
                        end do
                        !
                        !     case m = 1 through nlat-1
                        !
                        if (mmax < 2) exit case_block

                        do mp1=2, mmax
                            m = mp1-1
                            mb = m*(nlat-1)-(m*(m-1))/2
                            mp2 = mp1+1
                            if (mp1 <= ndo1) then
                                do k=1, nt
                                    do np1=mp1, ndo1, 2
                                        mn = mb+np1
                                        do i=1, imm1
                                            ve(i, 2*mp1-2, k) = ve(i, 2*mp1-2, k)-ci(mp1, np1, k)*wb(i, mn)
                                            ve(i, 2*mp1-1, k) = ve(i, 2*mp1-1, k)+cr(mp1, np1, k)*wb(i, mn)
                                            wo(i, 2*mp1-2, k) = wo(i, 2*mp1-2, k)-cr(mp1, np1, k)*vb(i, mn)
                                            wo(i, 2*mp1-1, k) = wo(i, 2*mp1-1, k)-ci(mp1, np1, k)*vb(i, mn)
                                        end do
                                        if (mlat /= 0) then
                                            ve(imid, 2*mp1-2, k) = ve(imid, 2*mp1-2, k) &
                                                -ci(mp1, np1, k)*wb(imid, mn)
                                            ve(imid, 2*mp1-1, k) = ve(imid, 2*mp1-1, k) &
                                                +cr(mp1, np1, k)*wb(imid, mn)
                                        end if
                                    end do
                                end do
                            end if

                            if (mp2 <= ndo2) then
                                do k=1, nt
                                    do np1=mp2, ndo2, 2
                                        mn = mb+np1
                                        do i=1, imm1
                                            ve(i, 2*mp1-2, k) = ve(i, 2*mp1-2, k)+br(mp1, np1, k)*vb(i, mn)
                                            ve(i, 2*mp1-1, k) = ve(i, 2*mp1-1, k)+bi(mp1, np1, k)*vb(i, mn)
                                            wo(i, 2*mp1-2, k) = wo(i, 2*mp1-2, k)-bi(mp1, np1, k)*wb(i, mn)
                                            wo(i, 2*mp1-1, k) = wo(i, 2*mp1-1, k)+br(mp1, np1, k)*wb(i, mn)
                                        end do
                                        if (mlat /= 0) then
                                            ve(imid, 2*mp1-2, k) = ve(imid, 2*mp1-2, k) &
                                                +br(mp1, np1, k)*vb(imid, mn)
                                            ve(imid, 2*mp1-1, k) = ve(imid, 2*mp1-1, k) &
                                                +bi(mp1, np1, k)*vb(imid, mn)
                                        end if
                                    end do
                                end do
                            end if
                        end do
                    case (4)
                        !
                        !     case ityp=4   v even,  w odd, and both cr and ci equal zero
                        !
                        !     case m = 0
                        !
                        do k=1, nt
                            do np1=2, ndo2, 2
                                do i=1, imid
                                    ve(i, 1, k)=ve(i, 1, k)+br(1, np1, k)*vb(i, np1)
                                end do
                            end do
                        end do
                        !
                        !     case m = 1 through nlat-1
                        !
                        if (mmax < 2) exit case_block

                        do mp1=2, mmax
                            m = mp1-1
                            mb = m*(nlat-1)-(m*(m-1))/2
                            mp2 = mp1+1
                            if (mp2 <= ndo2) then
                                do k=1, nt
                                    do np1=mp2, ndo2, 2
                                        mn = mb+np1
                                        do i=1, imm1
                                            ve(i, 2*mp1-2, k) = ve(i, 2*mp1-2, k)+br(mp1, np1, k)*vb(i, mn)
                                            ve(i, 2*mp1-1, k) = ve(i, 2*mp1-1, k)+bi(mp1, np1, k)*vb(i, mn)
                                            wo(i, 2*mp1-2, k) = wo(i, 2*mp1-2, k)-bi(mp1, np1, k)*wb(i, mn)
                                            wo(i, 2*mp1-1, k) = wo(i, 2*mp1-1, k)+br(mp1, np1, k)*wb(i, mn)
                                        end do
                                        if (mlat /= 0) then
                                            ve(imid, 2*mp1-2, k) = ve(imid, 2*mp1-2, k) &
                                                +br(mp1, np1, k)*vb(imid, mn)
                                            ve(imid, 2*mp1-1, k) = ve(imid, 2*mp1-1, k) &
                                                +bi(mp1, np1, k)*vb(imid, mn)
                                        end if
                                    end do
                                end do
                            end if
                        end do
                    case (5)
                        !
                        !     case ityp=5   v even,  w odd,     br and bi equal zero
                        !
                        !     case m = 0
                        !
                        do k=1, nt
                            do np1=3, ndo1, 2
                                do i=1, imm1
                                    wo(i, 1, k)=wo(i, 1, k)-cr(1, np1, k)*vb(i, np1)
                                end do
                            end do
                        end do
                        !
                        !     case m = 1 through nlat-1
                        !
                        if (mmax < 2) exit case_block

                        do mp1=2, mmax
                            m = mp1-1
                            mb = m*(nlat-1)-(m*(m-1))/2
                            mp2 = mp1+1

                            if (mp1 <= ndo1) then
                                do k=1, nt
                                    do np1=mp1, ndo1, 2
                                        mn = mb+np1
                                        do i=1, imm1
                                            ve(i, 2*mp1-2, k) = ve(i, 2*mp1-2, k)-ci(mp1, np1, k)*wb(i, mn)
                                            ve(i, 2*mp1-1, k) = ve(i, 2*mp1-1, k)+cr(mp1, np1, k)*wb(i, mn)
                                            wo(i, 2*mp1-2, k) = wo(i, 2*mp1-2, k)-cr(mp1, np1, k)*vb(i, mn)
                                            wo(i, 2*mp1-1, k) = wo(i, 2*mp1-1, k)-ci(mp1, np1, k)*vb(i, mn)
                                        end do
                                        if (mlat /= 0) then
                                            ve(imid, 2*mp1-2, k) = ve(imid, 2*mp1-2, k) &
                                                -ci(mp1, np1, k)*wb(imid, mn)
                                            ve(imid, 2*mp1-1, k) = ve(imid, 2*mp1-1, k) &
                                                +cr(mp1, np1, k)*wb(imid, mn)
                                        end if
                                    end do
                                end do
                            end if
                        end do
                    case (6)
                        !
                        !     case ityp=6   v odd  ,  w even
                        !
                        !     case m = 0
                        !
                        do k=1, nt
                            do np1=2, ndo2, 2
                                do i=1, imid
                                    we(i, 1, k)=we(i, 1, k)-cr(1, np1, k)*vb(i, np1)
                                end do
                            end do
                        end do

                        do k=1, nt
                            do np1=3, ndo1, 2
                                do i=1, imm1
                                    vo(i, 1, k)=vo(i, 1, k)+br(1, np1, k)*vb(i, np1)
                                end do
                            end do
                        end do
                        !
                        !     case m = 1 through nlat-1
                        !
                        if (mmax < 2) exit case_block

                        do mp1=2, mmax
                            m = mp1-1
                            mb = m*(nlat-1)-(m*(m-1))/2
                            mp2 = mp1+1
                            if (mp1 <= ndo1) then
                                do k=1, nt
                                    do np1=mp1, ndo1, 2
                                        mn = mb+np1
                                        do i=1, imm1
                                            vo(i, 2*mp1-2, k) = vo(i, 2*mp1-2, k)+br(mp1, np1, k)*vb(i, mn)
                                            vo(i, 2*mp1-1, k) = vo(i, 2*mp1-1, k)+bi(mp1, np1, k)*vb(i, mn)
                                            we(i, 2*mp1-2, k) = we(i, 2*mp1-2, k)-bi(mp1, np1, k)*wb(i, mn)
                                            we(i, 2*mp1-1, k) = we(i, 2*mp1-1, k)+br(mp1, np1, k)*wb(i, mn)
                                        end do
                                        if (mlat /= 0) then
                                            we(imid, 2*mp1-2, k) = we(imid, 2*mp1-2, k) &
                                                -bi(mp1, np1, k)*wb(imid, mn)
                                            we(imid, 2*mp1-1, k) = we(imid, 2*mp1-1, k) &
                                                +br(mp1, np1, k)*wb(imid, mn)
                                        end if
                                    end do
                                end do
                            end if

                            if (mp2 <= ndo2) then
                                do k=1, nt
                                    do np1=mp2, ndo2, 2
                                        mn = mb+np1
                                        do i=1, imm1
                                            vo(i, 2*mp1-2, k) = vo(i, 2*mp1-2, k)-ci(mp1, np1, k)*wb(i, mn)
                                            vo(i, 2*mp1-1, k) = vo(i, 2*mp1-1, k)+cr(mp1, np1, k)*wb(i, mn)
                                            we(i, 2*mp1-2, k) = we(i, 2*mp1-2, k)-cr(mp1, np1, k)*vb(i, mn)
                                            we(i, 2*mp1-1, k) = we(i, 2*mp1-1, k)-ci(mp1, np1, k)*vb(i, mn)
                                        end do
                                        if (mlat /= 0) then
                                            we(imid, 2*mp1-2, k) = we(imid, 2*mp1-2, k) &
                                                -cr(mp1, np1, k)*vb(imid, mn)
                                            we(imid, 2*mp1-1, k) = we(imid, 2*mp1-1, k) &
                                                -ci(mp1, np1, k)*vb(imid, mn)
                                        end if
                                    end do
                                end do
                            end if
                        end do
                    case (7)
                        !
                        !     case ityp=7   v odd, w even   cr and ci equal zero
                        !
                        !     case m = 0
                        !
                        do k=1, nt
                            do np1=3, ndo1, 2
                                do i=1, imm1
                                    vo(i, 1, k)=vo(i, 1, k)+br(1, np1, k)*vb(i, np1)
                                end do
                            end do
                        end do
                        !
                        !     case m = 1 through nlat-1
                        !
                        if (mmax < 2) exit case_block

                        do mp1=2, mmax
                            m = mp1-1
                            mb = m*(nlat-1)-(m*(m-1))/2
                            mp2 = mp1+1
                            if (mp1 <= ndo1) then
                                do k=1, nt
                                    do np1=mp1, ndo1, 2
                                        mn = mb+np1
                                        do i=1, imm1
                                            vo(i, 2*mp1-2, k) = vo(i, 2*mp1-2, k)+br(mp1, np1, k)*vb(i, mn)
                                            vo(i, 2*mp1-1, k) = vo(i, 2*mp1-1, k)+bi(mp1, np1, k)*vb(i, mn)
                                            we(i, 2*mp1-2, k) = we(i, 2*mp1-2, k)-bi(mp1, np1, k)*wb(i, mn)
                                            we(i, 2*mp1-1, k) = we(i, 2*mp1-1, k)+br(mp1, np1, k)*wb(i, mn)
                                        end do
                                        if (mlat /= 0) then
                                            we(imid, 2*mp1-2, k) = we(imid, 2*mp1-2, k) &
                                                -bi(mp1, np1, k)*wb(imid, mn)
                                            we(imid, 2*mp1-1, k) = we(imid, 2*mp1-1, k) &
                                                +br(mp1, np1, k)*wb(imid, mn)
                                        end if
                                    end do
                                end do
                            end if
                        end do
                    case (8)
                        !
                        !     case ityp=8   v odd,  w even   br and bi equal zero
                        !
                        !     case m = 0
                        !
                        do k=1, nt
                            do np1=2, ndo2, 2
                                do i=1, imid
                                    we(i, 1, k)=we(i, 1, k)-cr(1, np1, k)*vb(i, np1)
                                end do
                            end do
                        end do
                        !
                        !     case m = 1 through nlat-1
                        !
                        if (mmax < 2) exit case_block

                        do mp1=2, mmax
                            m = mp1-1
                            mb = m*(nlat-1)-(m*(m-1))/2
                            mp2 = mp1+1
                            if (mp2 <= ndo2) then
                                do k=1, nt
                                    do np1=mp2, ndo2, 2
                                        mn = mb+np1
                                        do i=1, imm1
                                            vo(i, 2*mp1-2, k) = vo(i, 2*mp1-2, k)-ci(mp1, np1, k)*wb(i, mn)
                                            vo(i, 2*mp1-1, k) = vo(i, 2*mp1-1, k)+cr(mp1, np1, k)*wb(i, mn)
                                            we(i, 2*mp1-2, k) = we(i, 2*mp1-2, k)-cr(mp1, np1, k)*vb(i, mn)
                                            we(i, 2*mp1-1, k) = we(i, 2*mp1-1, k)-ci(mp1, np1, k)*vb(i, mn)
                                        end do
                                        if (mlat /= 0) then
                                            we(imid, 2*mp1-2, k) = we(imid, 2*mp1-2, k) &
                                                -cr(mp1, np1, k)*vb(imid, mn)
                                            we(imid, 2*mp1-1, k) = we(imid, 2*mp1-1, k) &
                                                -ci(mp1, np1, k)*vb(imid, mn)
                                        end if
                                    end do
                                end do
                            end if
                        end do
                end select
            end block case_block

            do k=1, nt
                call hfft%backward(idv, nlon, ve(1, 1, k), idv, wrfft, work)
                call hfft%backward(idv, nlon, we(1, 1, k), idv, wrfft, work)
            end do

            if (ityp <= 2) then
                do k=1, nt
                    do j=1, nlon
                        do i=1, imm1
                            v(i, j, k) = 0.5_wp*(ve(i, j, k)+vo(i, j, k))
                            w(i, j, k) = 0.5_wp*(we(i, j, k)+wo(i, j, k))
                            v(nlp1-i, j, k) = 0.5_wp*(ve(i, j, k)-vo(i, j, k))
                            w(nlp1-i, j, k) = 0.5_wp*(we(i, j, k)-wo(i, j, k))
                        end do
                    end do
                end do
            else
                do k=1, nt
                    do j=1, nlon
                        v(1: imm1, j, k) = 0.5_wp*ve(1: imm1, j, k)
                        w(1: imm1, j, k) = 0.5_wp*we(1: imm1, j, k)
                    end do
                end do
            end if

            if (mlat /= 0) then
                do k=1, nt
                    v(imid, 1: nlon, k) = 0.5_wp*ve(imid, 1: nlon, k)
                    w(imid, 1: nlon, k) = 0.5_wp*we(imid, 1: nlon, k)
                end do
            end if

        end subroutine vhses1

    end subroutine vhses



    subroutine vhsesi(nlat, nlon, wvhses, lvhses, work, lwork, dwork, &
        ldwork, ierror)
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        integer (ip), intent (in)  :: nlat
        integer (ip), intent (in)  :: nlon
        real (wp),    intent (out) :: wvhses(lvhses)
        integer (ip), intent (in)  :: lvhses
        real (wp),    intent (out) :: work(lwork)
        integer (ip), intent (in)  :: lwork
        real (wp),    intent (out) :: dwork(ldwork)
        integer (ip), intent (in)  :: ldwork
        integer (ip), intent (out) :: ierror
        !----------------------------------------------------------------------
        ! Dictionary: local variables
        !----------------------------------------------------------------------
        integer (ip)    :: imid, labc, lzimn, mmax
        integer (ip)    :: workspace_indices(4)
        type (HFFTpack) :: hfft
        !----------------------------------------------------------------------

        mmax = min(nlat, (nlon+1)/2)
        imid = (nlat+1)/2
        lzimn = (imid*mmax*(2*nlat-mmax+1))/2
        labc = 3*(max(mmax-2, 0)*(2*nlat-mmax-1))/2

        !
        !==> Check validity of input arguments
        !
        if (nlat < 3) then
            ierror = 1
            return
        else if (nlon < 1) then
            ierror = 2
            return
        else if (lvhses < 2*lzimn+nlon+15) then
            ierror = 3
            return
        else if (lwork < 5*nlat*imid+labc) then
            ierror = 4
            return
        else if (ldwork < 2*(nlat+1)) then
            ierror = 5
            return
        else
            ierror = 0
        end if

        !
        !==> Set workspace indices
        !
        workspace_indices = get_workspace_indices(nlat, imid, mmax, lzimn)

        associate( &
            iw1 => workspace_indices(1), &
            idz => workspace_indices(2), &
            jw1 => workspace_indices(3), &
            jw2 => workspace_indices(4) &
            )

            call ves1(nlat, nlon, imid, wvhses, wvhses(jw1), idz, work, work(iw1), dwork)

            call hfft%initialize(nlon, wvhses(jw2))

        end associate

    contains


        pure function get_workspace_indices(nlat, imid, mmax, lzimn) result (return_value)
            !----------------------------------------------------------------------
            ! Dictionary: calling arguments
            !----------------------------------------------------------------------
            integer (ip), intent (in)  :: nlat
            integer (ip), intent (in)  :: imid
            integer (ip), intent (in)  :: mmax
            integer (ip), intent (in)  :: lzimn
            integer (ip)               :: return_value(4)
            !----------------------------------------------------------------------

            associate( i => return_value )

                i(1) = 3*nlat*imid+1
                i(2) = (mmax*(2*nlat-mmax+1))/2
                i(3) = lzimn+1
                i(4) = 2*lzimn+1

            end associate

        end function get_workspace_indices


        subroutine ves1(nlat, nlon, imid, vb, wb, idz, vin, wzvin, dwork)
            !----------------------------------------------------------------------
            ! Dictionary: calling arguments
            !----------------------------------------------------------------------
            integer (ip), intent (in)  :: nlat
            integer (ip), intent (in)  :: nlon
            integer (ip), intent (in)  :: imid
            real (wp),    intent (out) :: vb(imid, *)
            real (wp),    intent (out) :: wb(imid, *)
            integer (ip), intent (in)  :: idz
            real (wp),    intent (out) :: vin(imid, nlat, 3)
            real (wp),    intent (out) :: wzvin(*)
            real (wp),    intent (out) :: dwork(*)
            !----------------------------------------------------------------------
            ! Dictionary: local variables
            !----------------------------------------------------------------------
            integer (ip)         :: i3, m, mmax, mn, mp1, np1
            type (SpherepackAux) :: sphere_aux
            !----------------------------------------------------------------------

            mmax = min(nlat, (nlon+1)/2)

            call sphere_aux%vbinit(nlat, nlon, wzvin, dwork)

            do mp1=1, mmax
                m = mp1-1
                call sphere_aux%vbin(0, nlat, nlon, m, vin, i3, wzvin)
                do np1=mp1, nlat
                    mn = m*(nlat-1)-(m*(m-1))/2+np1
                    vb(1: imid, mn) = vin(1: imid, np1, i3)
                end do
            end do

            call sphere_aux%wbinit(nlat, nlon, wzvin, dwork)

            do mp1=1, mmax
                m = mp1-1
                call sphere_aux%wbin(0, nlat, nlon, m, vin, i3, wzvin)
                do np1=mp1, nlat
                    mn = m*(nlat-1)-(m*(m-1))/2+np1
                    wb(1: imid, mn) = vin(1: imid, np1, i3)
                end do
            end do

        end subroutine ves1

    end subroutine vhsesi

end module module_vhses
